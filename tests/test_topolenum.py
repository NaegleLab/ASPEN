import unittest
import os, os.path
import time
import cPickle as pickle
from cStringIO import StringIO
from multiprocessing import Condition,Pipe
from mock import patch,call,mock_open,Mock,PropertyMock
from aspen import topolenum as te


@patch('tempfile.NamedTemporaryFile')
class TestFIFOfileBaseClassTMPFILEclass(unittest.TestCase):
  
  def setUp(self):
    reload(te)
    self.TMPFILE = te.FIFOfile.TMPFILE
    self.TMPFILE.init_class(mode='b',wbuf=0,rbuf=0,suffix='',delete=True,
                            dir='/dummy/path',check_delay=5)
  
  def test_instance_counting_and_file_creation_and_naming(self,patched_NTF):
    
    tmpfile_instance1 = self.TMPFILE()
    self.assertEqual(tmpfile_instance1.instcount,self.TMPFILE.instcount)
    self.assertEqual(self.TMPFILE.instcount,1)
    patched_NTF.assert_called_with('rb',0,suffix='',prefix='FIFOfile001_',
                                   dir='/dummy/path',delete=True)
    self.TMPFILE()
    self.assertEqual(self.TMPFILE.instcount,2)
    patched_NTF.assert_called_with('rb',0,suffix='',prefix='FIFOfile002_',
                                   dir='/dummy/path',delete=True)
  
  def test_spooling_and_file_spool_interaction(self,patched_NTF):
    self.assertFalse(hasattr(self.TMPFILE,'file_spool'))
    self.TMPFILE.start_spooling()
    self.assertTrue(hasattr(self.TMPFILE,'file_spool'))
    self.assertEqual(len(self.TMPFILE.file_spool),0)
    tmpfile_instance = self.TMPFILE.spool()
    self.assertEqual(len(self.TMPFILE.file_spool),1)
    self.assertIs(tmpfile_instance,self.TMPFILE.pop_from_spool())
    self.assertEqual(len(self.TMPFILE.file_spool),0)
    self.assertIs(None,self.TMPFILE.pop_from_spool())
  
  def test_writing_handle(self,patched_NTF):
    tmpfile_instance = self.TMPFILE()
    self.assertFalse(hasattr(tmpfile_instance,'wh'))
    with patch('__builtin__.open',mock_open(),create=True) as patched_open:
      tmpfile_instance.open()
      patched_open.assert_called_once_with(patched_NTF.return_value.name,'wb',0)
      self.assertTrue(hasattr(tmpfile_instance,'wh'))
      self.assertIs(tmpfile_instance.wh,patched_open.return_value)
    tmpfile_instance.close()
    tmpfile_instance.wh.close.assert_called_once_with()
  
  @patch('os.path.getsize',side_effect=[10.0,20.0])
  def test_size_checking(self,patched_getsize,patched_NTF):
    tmpfile_instance = self.TMPFILE()
    self.assertSequenceEqual([tmpfile_instance.size for i in xrange(5)],
                             [0.0,0.0,0.0,0.0,0.0])
    patched_getsize.assert_not_called()
    self.assertEqual(tmpfile_instance.size,10.0)
    patched_getsize.assert_called_once_with(patched_NTF.return_value.name)
    self.assertSequenceEqual([tmpfile_instance.size for i in xrange(5)],
                             [10.0,10.0,10.0,10.0,10.0])
    patched_getsize.assert_called_once_with(patched_NTF.return_value.name)
    self.assertEqual(tmpfile_instance.size,20.0)
    self.assertSequenceEqual(patched_getsize.call_args_list,
                             [call(patched_NTF.return_value.name),
                              call(patched_NTF.return_value.name)],list)
  
  @patch('os.path.exists',return_value=False)
  def test_file_discarding(self,patched_exists,patched_NTF):
    tmpfile_instance = self.TMPFILE()
    self.assertIs(tmpfile_instance.rh,patched_NTF.return_value)
    tmpfile_instance.discard()
    tmpfile_instance.rh.close.assert_called_once_with()


@patch('os.path.exists',return_value=False)
@patch('os.path.realpath',return_value='/dummy/path')
@patch('aspen.topolenum.TemporaryDirectory')
@patch('tempfile.NamedTemporaryFile',**{'return_value.name':'dummy_temp_file'})
class TestFIFOfileBaseClass(unittest.TestCase):
  
  def setUp(self):
    reload(te)
  
  @patch('aspen.topolenum.FIFOfile.TMPFILE.init_class')
  def test_tmpdir_creation_and_TMPFILE_class_init(self,patched_TMPFILE_initcls,
                                                  patched_NTF,patched_TmpDir,
                                                  patched_realpath,
                                                  patched_exists):
    fifo_obj = te.FIFOfile(top_path='/dummy/top/path/arg',suffix='dummy_suffix')
    patched_TmpDir.assert_called_once_with(dir='/dummy/top/path/arg',
                                           prefix='FIFOworkspace_',
                                           suffix='dummy_suffix')
    patched_TMPFILE_initcls.assert_called_once_with('b',0,0,'dummy_suffix',True,
                                                    '/dummy/path',100)
  
  def test_starting_fifo_OUT_end(self,patched_NTF,*args):
    fifo_obj = te.FIFOfile()
    self.assertFalse(hasattr(te.FIFOfile.TMPFILE,'file_spool'))
    patched_NTF.assert_not_called()
    self.assertFalse(hasattr(fifo_obj,'current_reading_file'))
    
    fifo_obj.start_OUT_end()
    patched_NTF.assert_called_once_with('rb',0,suffix='',prefix='FIFOfile001_',
                                        dir='/dummy/path',delete=True)
    self.assertSequenceEqual(te.FIFOfile.TMPFILE.file_spool,[],seq_type=list)
    self.assertEqual(fifo_obj.current_reading_file.name,'dummy_temp_file')
  
  def test_starting_fifo_IN_end(self,*args):
    fifo_obj = te.FIFOfile()
    self.assertFalse(hasattr(fifo_obj,'current_writing_file'))
    
    fifo_obj.start_OUT_end()
    with patch('__builtin__.open',mock_open(),create=True) as patched_open:
      fifo_obj.start_IN_end()
      self.assertIs(fifo_obj.current_writing_file,fifo_obj.current_reading_file)
      patched_open.assert_called_once_with('dummy_temp_file','wb',0)
      self.assertIs(fifo_obj.current_writing_file.wh,patched_open.return_value)
  
  @patch('os.path.getsize',side_effect=[500.0,1100.0])
  @patch('__builtin__.open')
  def test_wh_retrieval_and_rollover(self,patched_open,patched_getsize,
                                     patched_NTF,*args):
    def patched_NTF_side_effect(*args,**kwargs):
      return_val = Mock()
      return_val.name = kwargs['prefix']
      return return_val
    patched_NTF.side_effect = patched_NTF_side_effect
    
    def assert_no_rollover():
      self.assertFalse(fifo_obj.current_reading_file.wh.close.called)
      self.assertFalse(patched_NTF.called)
      self.assertFalse(patched_open.called)
      self.assertIs(fifo_obj.current_writing_file,
                    fifo_obj.current_reading_file)
      self.assertSequenceEqual(fifo_obj.TMPFILE.file_spool,[],seq_type=list)
    
    # Start up FIFO and clear out calls to mocks
    # 0.001*1024^3 ~= 1073.7, 500.0 < 1073.7 < 1100.0, perfect!
    fifo_obj = te.FIFOfile(max_file_size_GB=0.000001,size_check_delay=2)
    fifo_obj.start_OUT_end()
    fifo_obj.start_IN_end()
    patched_NTF.reset_mock()
    patched_open.reset_mock()
    
    # First and second calls should trigger nothing
    fifo_obj.wh
    fifo_obj.wh
    self.assertFalse(patched_getsize.called)
    assert_no_rollover()
    
    # Third call should trigger file size check ...
    fifo_obj.wh
    patched_getsize.assert_called_once_with('FIFOfile001_')
    patched_getsize.reset_mock()
    self.assertEqual(fifo_obj.current_writing_file._size,500.0)
    # ... but the file size should be insufficient to trigger a rollover
    assert_no_rollover()
    
    # Fourth and fifth calls should, again, trigger nothing
    fifo_obj.wh
    fifo_obj.wh
    self.assertFalse(patched_getsize.called)
    assert_no_rollover()
    
    # Finally, sixth call should trigger a second size check ...
    fifo_obj.wh
    patched_getsize.assert_called_once_with('FIFOfile001_')
    # ... and the returned file size should trigger a rollover to a new file, ...
    fifo_obj.current_reading_file.wh.close.assert_called_once_with()
    patched_NTF.assert_called_once_with('rb',0,prefix='FIFOfile002_',suffix='',
                                        dir='/dummy/path',delete=True)
    patched_open.assert_called_once_with('FIFOfile002_','wb',0)
    self.assertSequenceEqual(fifo_obj.TMPFILE.file_spool,
                             [fifo_obj.current_writing_file])
    # ... meaning the reading and writing files are now different ...
    self.assertIsNot(fifo_obj.current_writing_file,
                     fifo_obj.current_reading_file)
    self.assertEqual(fifo_obj.current_writing_file._size,0.0)
    self.assertEqual(fifo_obj.current_reading_file._size,1100.0)
  
  def test_rh_retrival_EOF_testing_and_discarding(self,patched_NTF,*args):
    patched_NTF.instances = []
    def patched_NTF_side_effect(*args,**kwargs):
      return_val = Mock(**{'tell.side_effect':[123,124,200,200,200]})
      return_val.name = kwargs['prefix']
      patched_NTF.instances.append(return_val)
      return return_val
    patched_NTF.side_effect = patched_NTF_side_effect
    
    # Start up FIFO
    fifo_obj = te.FIFOfile()
    fifo_obj.start_OUT_end()
    with patch('__builtin__.open',mock_open(),create=True):
      fifo_obj.start_IN_end()
    
    # On first retrieval attempt rh has something to read
    self.assertIs(fifo_obj.current_reading_file.rh,fifo_obj.rh)
    self.assertIs(patched_NTF.instances[0],fifo_obj.current_reading_file.rh)
    self.assertSequenceEqual(fifo_obj.current_reading_file.rh.tell.call_args_list,
                             [call(),call()])
    fifo_obj.current_reading_file.rh.read.assert_called_once_with(1)
    fifo_obj.current_reading_file.rh.seek.assert_called_once_with(123)
    
    # Spool second temp file to which rh can rollover
    fifo_obj.current_writing_file = fifo_obj.TMPFILE.spool()
    
    # On second retrieval attempt rh is at EOF - cursor not advanced after read()
    # At this point rh should rollover to next file
    self.assertIsNot(fifo_obj.current_reading_file.rh,fifo_obj.rh)
    self.assertIs(patched_NTF.instances[1],fifo_obj.current_reading_file.rh)
    # After rollover to new file cursor is repositioned to first value seek()
    # seek() returns when called on that file's handle
    fifo_obj.current_reading_file.rh.seek.assert_called_with(123)
    fifo_obj.current_reading_file.rh.seek.reset_mock()
    patched_NTF.instances[0].close.assert_called_once_with()
    
    # On third retrieval attempt the rh of the new file has something to read
    self.assertIs(fifo_obj.current_reading_file.rh,fifo_obj.rh)
    self.assertIs(patched_NTF.instances[1],fifo_obj.current_reading_file.rh)
    # Value 123 was already used up on previous retrieval of rh, so calls to
    # seek() before and after a call to read() return 124 and 200 - still
    # indicative of having something to read. Cursor is then returned to 124.
    fifo_obj.current_reading_file.rh.seek.assert_called_once_with(124)
    fifo_obj.current_reading_file.rh.seek.reset_mock()
    
    # On fourth retrieval attempt rh has nothing to read again, but also nothing
    # to rollover to, so despite being at EOF the same handle is returned
    self.assertIs(fifo_obj.current_reading_file.rh,fifo_obj.rh)
    self.assertFalse(fifo_obj.current_reading_file.rh.close.called)
    # Calls to seek() returned 200 and 200, indicating EOF. A seek(200) call is
    # still made to make sure handle is readable if new data is written.
    fifo_obj.current_reading_file.rh.seek.assert_called_once_with(200)
  
  def test_FIFO_pop(self,patched_NTF,patched_TmpDir,patched_realpath,
                    patched_exists):
    with patch.object(te.FIFOfile,'rh',PropertyMock(side_effect=\
                        [StringIO(pickle.dumps('data1',protocol=2)),
                         StringIO(''),
                         StringIO(pickle.dumps('data2',protocol=2))]))\
                                                                as patched_rh:
      fifo_obj = te.FIFOfile()
      fifo_obj.start_OUT_end()
      self.assertEqual(fifo_obj.pop(),'data1')
      self.assertEqual(fifo_obj.pop(),None)
      self.assertEqual(fifo_obj.pop(),'data2')
  
  @patch('__builtin__.open')
  def test_FIFO_close(self,patched_open,patched_NTF,patched_TmpDir,
                      patched_realpath,patched_exists):
    # Simple case: one temp file, spool is empty
    fifo_obj = te.FIFOfile()
    fifo_obj.start_OUT_end()
    fifo_obj.start_IN_end()
    
    self.assertIs(fifo_obj.tmpdir_obj,patched_TmpDir.return_value)
    self.assertIs(fifo_obj.current_reading_file,fifo_obj.current_writing_file)
    self.assertIs(fifo_obj.current_reading_file.rh,patched_NTF.return_value)
    self.assertIs(fifo_obj.current_writing_file.wh,patched_open.return_value)
    fifo_obj.close()
    patched_open.return_value.close.assert_called_once_with()
    patched_NTF.return_value.close.assert_called_once_with()
    patched_exists.assert_called_once_with('dummy_temp_file')
    fifo_obj.tmpdir_obj.__exit__.assert_called_once_with(None,None,None)
    
    def prepare_return_val(name,append_here):
      return_val = Mock()
      return_val.name = name
      return_val.closed = False
      append_here.append(return_val)
      return return_val
    
    patched_NTF.instances = []
    def patched_NTF_side_effect(*args,**kwargs):
      return prepare_return_val(kwargs['prefix'],patched_NTF.instances)
    patched_NTF.side_effect = patched_NTF_side_effect
    
    patched_open.instances = []
    def patched_open_side_effect(*args,**kwargs):
      return prepare_return_val(args[0],patched_open.instances)
    patched_open.side_effect = patched_open_side_effect
    
    # Alternative: reading and writing files are different, spool not empty
    fifo_obj = te.FIFOfile()
    fifo_obj.start_OUT_end()
    fifo_obj.start_IN_end()
    fifo_obj.current_writing_file.close()
    for i in xrange(2):
      fifo_obj.current_writing_file.wh.closed = True
      fifo_obj.current_writing_file = fifo_obj.TMPFILE.spool()
      fifo_obj.current_writing_file.open()
    
    self.assertIs(fifo_obj.tmpdir_obj,patched_TmpDir.return_value)
    patched_TmpDir.return_value.reset_mock()
    self.assertEqual(len(patched_NTF.instances),3)
    self.assertEqual(len(patched_open.instances),3)
    self.assertIsNot(fifo_obj.current_reading_file,
                     fifo_obj.current_writing_file)
    self.assertIs(fifo_obj.current_reading_file.rh,patched_NTF.instances[0])
    self.assertIs(fifo_obj.current_writing_file.wh,patched_open.instances[2])
    for i in xrange(3):
      if i < 2:
        self.assertIs(fifo_obj.TMPFILE.file_spool[i].rh,patched_NTF.instances[i+1])
        self.assertIs(fifo_obj.TMPFILE.file_spool[i].wh,patched_open.instances[i+1])
      self.assertEqual(patched_NTF.instances[i].name,
                       patched_open.instances[i].name)
    fifo_obj.close()
    for i in xrange(3):
      patched_open.instances[i].close.assert_called_once_with()
      patched_NTF.instances[i].close.assert_called_once_with()
    patched_exists.assert_has_calls([call('FIFOfile001_'),call('FIFOfile002_'),
                                     call('FIFOfile003_')])
    fifo_obj.tmpdir_obj.__exit__.assert_called_once_with(None,None,None)


class TestFIFOfile_and_TMPFILE_and_system_integration(unittest.TestCase):
  
  def setUp(self):
    reload(te)
    self.fifo_obj = te.FIFOfile(top_path=None, # Let os decide where to put tempdir
                                # 1.7695128917694092e-08 GB * 1024^4 B/GB = 19.0 B
                                max_file_size_GB=1.7695128917694092e-08,
                                size_check_delay=0) # Check file size every time
    self.fifo_obj.start_OUT_end()
    self.fifo_obj.start_IN_end()
  
  def test_rollover_and_file_discarding(self):
    # After init there is one temp file and it is empty
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.assertEqual(os.path.getsize(self.fifo_obj.current_writing_file.name),0)
    
    # Pickling a single int increases file size by 5 bytes
    # After pushing four ints the file size is 20, which is > max_size == 19,
    # but prior to previous push file size was still < max_size, so wh has not
    # yet rolled over to a new file
    self.fifo_obj.push(1);self.fifo_obj.push(1);self.fifo_obj.push(1);self.fifo_obj.push(1)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.assertEqual(os.path.getsize(self.fifo_obj.current_writing_file.name),20)
    
    # When wh is retrieved for next push it will rollover to a new file
    prev_writing_file_name = self.fifo_obj.current_writing_file.name
    self.fifo_obj.push(1)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),2)
    self.assertNotEqual(self.fifo_obj.current_writing_file.name,
                        prev_writing_file_name)
    self.assertEqual(os.path.getsize(prev_writing_file_name),20)
    self.assertEqual(os.path.getsize(self.fifo_obj.current_writing_file.name),5)
    
    # The last of four pushes in a row should, again, cause a file rollover
    two_ago_writing_file_name = prev_writing_file_name
    prev_writing_file_name = self.fifo_obj.current_writing_file.name
    self.fifo_obj.push(1)
    self.fifo_obj.push(1)
    self.fifo_obj.push(1)
    self.fifo_obj.push(1)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),3)
    self.assertNotEqual(self.fifo_obj.current_writing_file.name,
                        prev_writing_file_name)
    self.assertNotEqual(self.fifo_obj.current_writing_file.name,
                        two_ago_writing_file_name)
    self.assertEqual(os.path.getsize(two_ago_writing_file_name),20)
    self.assertEqual(os.path.getsize(prev_writing_file_name),20)
    self.assertEqual(os.path.getsize(self.fifo_obj.current_writing_file.name),5)
    
    # Since nothing has been read from the FIFO, rh is still on the first file
    # After four pops the rh cursor is at EOF of first file, but it has not yet
    # detected this fact and the file has not yet been discarded
    self.assertEqual(self.fifo_obj.current_reading_file.name,
                     two_ago_writing_file_name)
    self.fifo_obj.pop()
    self.fifo_obj.pop()
    self.fifo_obj.pop()
    self.fifo_obj.pop()
    self.assertEqual(self.fifo_obj.current_reading_file.name,
                     two_ago_writing_file_name)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),3)
    
    # After next pop rh rolls over to next file in the spool and the previous
    # file is discarded (deleted from disk)
    self.fifo_obj.pop();
    self.assertEqual(self.fifo_obj.current_reading_file.name,
                     prev_writing_file_name)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),2)
    self.assertNotIn(two_ago_writing_file_name,
                     os.listdir(self.fifo_obj.tmpdir_obj.name))
    
    # The last of four more pops should cause another rollover and discarding
    # of second file, leaving only the current writing file on disk
    self.fifo_obj.pop()
    self.fifo_obj.pop()
    self.fifo_obj.pop()
    self.fifo_obj.pop()
    self.assertEqual(self.fifo_obj.current_reading_file.name,
                     self.fifo_obj.current_writing_file.name)
    self.assertIs(self.fifo_obj.current_reading_file,
                  self.fifo_obj.current_writing_file)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.assertNotIn(prev_writing_file_name,
                     os.listdir(self.fifo_obj.tmpdir_obj.name))
    self.assertSequenceEqual(os.listdir(self.fifo_obj.tmpdir_obj.name),
                             [os.path.basename(
                                    self.fifo_obj.current_writing_file.name)])
  
  def test_catching_up_to_wh_with_rollovers(self):
    # After four pushes and four pops rh cursor has caught up to wh, which has
    # not yet rolled over to a new file
    self.fifo_obj.push(1)
    self.fifo_obj.push(2)
    self.fifo_obj.push(3)
    self.fifo_obj.push(4)
    self.assertSequenceEqual([self.fifo_obj.pop() for i in xrange(4)],[1,2,3,4])
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    
    # Without further writing, subsequent pops return None, still only one file
    self.assertSequenceEqual([self.fifo_obj.pop() for i in xrange(4)],[None,
                                                                       None,
                                                                       None,
                                                                       None])
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    
    # After six more pushes the wh has rolled over twice, but rh is sill on
    # the first file
    self.fifo_obj.push(5);self.fifo_obj.push(6);self.fifo_obj.push(7)
    self.fifo_obj.push(8);self.fifo_obj.push(9);self.fifo_obj.push(10)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),3)
    
    # After six pops rh catches up to wh again, across two files
    self.assertSequenceEqual([self.fifo_obj.pop() for i in xrange(6)],
                             [5,6,7,8,9,10])
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.assertIs(self.fifo_obj.current_reading_file,
                  self.fifo_obj.current_writing_file)
    
    # Further pops again return Nones
    self.assertSequenceEqual([self.fifo_obj.pop() for i in xrange(4)],[None,
                                                                       None,
                                                                       None,
                                                                       None])
  
  def test_cleanup(self):
    # After pushing ten integers there are three temp files
    for i in xrange(10):
      self.fifo_obj.push(i+1)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),3)
    
    # After five pops there are two temp files remaining
    self.assertSequenceEqual([self.fifo_obj.pop() for i in xrange(5)],[1,2,3,4,5])
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),2)
    
    # After closing FIFO remaining files and the temp dir are gone from disk
    self.fifo_obj.close()
    self.assertFalse(os.path.exists(self.fifo_obj.current_writing_file.name))
    self.assertFalse(os.path.exists(self.fifo_obj.current_reading_file.name))
    self.assertFalse(os.path.exists(self.fifo_obj.tmpdir_obj.name))
  
  def tearDown(self):
    self.fifo_obj.close()


class TestSharedFIFOfileClassTMPFILEClass(unittest.TestCase):
   
  def setUp(self):
    reload(te)
    self.TMPFILE = te.SharedFIFOfile.TMPFILE
    self.TMPFILE.init_class(mode='b',wbuf=0,rbuf=0,suffix='',delete=True,
                            dir='/dummy/path',check_delay=0)
   
  @patch('os.path.getsize',return_value=0)
  def test_spool_call_and_writing_side_init(self,patched_getsize):
    self.TMPFILE.reading_side_conn.send('dummy_tempfile_name')
    tempfile_obj = self.TMPFILE.spool()
    self.assertIs(self.TMPFILE.reading_side_conn.recv(),None)
    patched_getsize.assert_called_once_with('dummy_tempfile_name')
    self.assertEqual(tempfile_obj.name,'dummy_tempfile_name')
    self.assertIs(tempfile_obj._size,patched_getsize.return_value)
    self.assertEqual(tempfile_obj.access_count_since_size_check,0)
  
  @patch('aspen.topolenum.FIFOfile.TMPFILE.__init__')
  def test_reading_side_init(self,patched_TMPFILE_init):
    tempfile_obj = self.TMPFILE()
    patched_TMPFILE_init.assert_called_once_with(tempfile_obj)


class TestSharedFIFOfileClassSpoolerThreadClass(unittest.TestCase):
  
  def setUp(self):
    reload(te)
    self.writing_end_conn,self.reading_end_conn = te.multiprocessing.Pipe()
    self.mock_spool_callable = Mock()
    self.tmpfile_name_mock_property = \
                              PropertyMock(return_value='dummy_tempfile_name')
    type(self.mock_spool_callable.return_value).name = \
                                                self.tmpfile_name_mock_property
    self.spooler = te.SharedFIFOfile.SpoolerThread(self.reading_end_conn,
                                                   self.mock_spool_callable,
                                                   interval_len=0.01)
  
  def test_SpoolerThread(self):
    self.assertEqual(self.spooler.name,'MainProcess--TmpfileSpoolerThread')
    self.assertEqual(len(te.threading.enumerate()),1)
    self.spooler.start()
    self.assertEqual(len(te.threading.enumerate()),2)
    self.assertFalse(self.writing_end_conn.poll())
    self.assertFalse(self.mock_spool_callable.called)
    self.assertFalse(self.tmpfile_name_mock_property.called)
    self.writing_end_conn.send(None)
    time.sleep(0.1) # Things don't happen instantaneously in the other thread
    self.mock_spool_callable.assert_called_once_with()
    self.tmpfile_name_mock_property.assert_called_once_with()
    self.assertTrue(self.writing_end_conn.poll())
    self.assertEqual(self.writing_end_conn.recv(),'dummy_tempfile_name')
  
  def tearDown(self):
    self.spooler.stop.set()
    self.spooler.join()


@patch('os.path.exists',return_value=False)
@patch('os.path.realpath',return_value='/dummy/path')
@patch('aspen.topolenum.TemporaryDirectory')
@patch('tempfile.NamedTemporaryFile',**{'return_value.name':'dummy_temp_file'})
class TestSharedFIFOfileClass(unittest.TestCase):
  
  def setUp(self):
    reload(te)
  
  def test_TMPFILE_class_init_w_conns(self,*args):
    self.assertFalse(hasattr(te.SharedFIFOfile.TMPFILE,'instcount'))
    self.assertFalse(hasattr(te.SharedFIFOfile.TMPFILE,'writing_side_conn'))
    self.assertFalse(hasattr(te.SharedFIFOfile.TMPFILE,'reading_side_conn'))
    
    fifo_obj = te.SharedFIFOfile()
    self.assertTrue(hasattr(te.SharedFIFOfile.TMPFILE,'instcount'))
    self.assertTrue(hasattr(te.SharedFIFOfile.TMPFILE,'writing_side_conn'))
    self.assertTrue(hasattr(te.SharedFIFOfile.TMPFILE,'reading_side_conn'))
    
    te.SharedFIFOfile.TMPFILE.writing_side_conn.send('test_packet1')
    self.assertEqual(te.SharedFIFOfile.TMPFILE.reading_side_conn.recv(),'test_packet1')
    te.SharedFIFOfile.TMPFILE.reading_side_conn.send('test_packet2')
    self.assertEqual(te.SharedFIFOfile.TMPFILE.writing_side_conn.recv(),'test_packet2')
  
  def test_starting_fifo_OUT_end(self,patched_NTF,*args):
    self.assertFalse(hasattr(te.SharedFIFOfile.TMPFILE,'file_spool'))
    self.assertEqual(len(te.threading.enumerate()),1)
    
    self.fifo_obj = te.SharedFIFOfile(interval_len=0.01)
    self.assertFalse(hasattr(self.fifo_obj,'side'))
    self.assertFalse(hasattr(self.fifo_obj,'current_reading_file'))
    self.assertFalse(hasattr(self.fifo_obj,'current_writing_file'))
    self.assertFalse(hasattr(self.fifo_obj,'spooler'))
    self.assertFalse(patched_NTF.called)
    
    self.fifo_obj.start_OUT_end()
    self.assertTrue(hasattr(self.fifo_obj,'side'))
    self.assertEqual(self.fifo_obj.side,'reading')
    
    self.assertEqual(patched_NTF.call_count,1)
    self.assertTrue(hasattr(self.fifo_obj,'current_reading_file'))
    self.assertFalse(hasattr(self.fifo_obj,'current_writing_file'))
    self.assertIs(self.fifo_obj.current_reading_file.rh,
                  patched_NTF.return_value)
    self.assertIs(self.fifo_obj.current_reading_file.rh.name,'dummy_temp_file')
    
    self.assertTrue(hasattr(te.SharedFIFOfile.TMPFILE,'file_spool'))
    self.assertSequenceEqual(te.SharedFIFOfile.TMPFILE.file_spool,[])
    
    self.assertTrue(te.SharedFIFOfile.TMPFILE.writing_side_conn.poll())
    self.assertEqual(te.SharedFIFOfile.TMPFILE.writing_side_conn.recv(),'dummy_temp_file')
    self.assertEqual(len(te.threading.enumerate()),2)
    self.assertTrue(hasattr(self.fifo_obj,'spooler'))
    self.assertTrue(self.fifo_obj.spooler.is_alive())
  
  @patch('os.path.getsize',return_value=100)
  def test_starting_fifo_IN_end(self,patched_getsize,*args):
    self.fifo_obj = te.SharedFIFOfile(interval_len=0.01)
    self.assertFalse(hasattr(self.fifo_obj,'side'))
    self.assertFalse(hasattr(self.fifo_obj,'current_reading_file'))
    self.assertFalse(hasattr(self.fifo_obj,'current_writing_file'))
    
    te.SharedFIFOfile.TMPFILE.reading_side_conn.send('dummy_temp_file')
    with patch('__builtin__.open',mock_open(),create=True) as patched_open:
      self.fifo_obj.start_IN_end()
    
    self.assertTrue(hasattr(self.fifo_obj,'side'))
    self.assertEqual(self.fifo_obj.side,'writing')
    
    self.assertFalse(hasattr(self.fifo_obj,'current_reading_file'))
    self.assertTrue(hasattr(self.fifo_obj,'current_writing_file'))
    self.assertEqual(self.fifo_obj.current_writing_file.name,'dummy_temp_file')
    patched_getsize.assert_called_once_with('dummy_temp_file')
    self.fifo_obj.current_writing_file._size == 100
    
    patched_open.assert_called_once_with('dummy_temp_file','wb',0)
    self.assertIs(self.fifo_obj.current_writing_file.wh,
                  patched_open.return_value)
  
  @patch('os.path.getsize',side_effect=[0,1100.0]*2)
  @patch('multiprocessing.Pipe')
  @patch('__builtin__.open')
  def test_wh_rollover_writing_end(self,patched_open,patched_Pipe,
                                   patched_getsize,patched_NTF,*args):
    def patched_Pipe_side_effect():
      writing_side_conn,reading_side_conn = Pipe()
      patched_Pipe.wrsconn = Mock(wraps=writing_side_conn)
      patched_Pipe.rdsconn = Mock(wraps=reading_side_conn)
      return patched_Pipe.wrsconn,patched_Pipe.rdsconn
    patched_Pipe.side_effect = patched_Pipe_side_effect
    
    def patched_NTF_side_effect(*args,**kwargs):
      return_val = Mock()
      return_val.name = kwargs['prefix']
      return return_val
    patched_NTF.side_effect = patched_NTF_side_effect
    
    class TesterProc(te.multiprocessing.Process):
      def __init__(self,fifo_obj,shutdown):
        te.multiprocessing.Process.__init__(self)
        self.fifo_obj = fifo_obj
        self.shutdown = shutdown
      
      def run(self):
        self.fifo_obj.start_OUT_end()
        while not self.shutdown.wait(0.01):
          pass
        self.fifo_obj.close()
    
    self.fifo_obj = te.SharedFIFOfile(max_file_size_GB=0.000001,
                                      size_check_delay=1,interval_len=0.01)
    self.shutdown_EV = te.multiprocessing.Event()
    self.other_proc = TesterProc(self.fifo_obj,self.shutdown_EV)
    self.other_proc.start()
    self.fifo_obj.start_IN_end()
    patched_Pipe.wrsconn.recv.reset_mock()
    patched_open.reset_mock()
    patched_getsize.reset_mock()
    
    self.fifo_obj.wh
    self.assertFalse(patched_getsize.called)
    self.assertFalse(patched_Pipe.wrsconn.send.called)
    self.assertFalse(patched_Pipe.wrsconn.recv.called)
    self.assertFalse(patched_open.return_value.close.called)
    self.assertFalse(patched_open.called)
    self.fifo_obj.wh
    patched_getsize.assert_has_calls([call('FIFOfile001_'),
                                      call('FIFOfile002_')])
    patched_Pipe.wrsconn.send.assert_called_once_with(None)
    patched_Pipe.wrsconn.recv.assert_called_once_with()
    patched_open.return_value.close.assert_called_once_with()
    patched_open.assert_called_once_with('FIFOfile002_','wb',0)
    
    self.shutdown_EV.set()
    self.fifo_obj.close()
    self.other_proc.join()
  
  @patch('os.path.getsize',side_effect=[0,1100.0]*2)
  @patch('multiprocessing.Pipe')
  @patch('__builtin__.open')
  def test_wh_rollover_reading_end(self,patched_open,patched_Pipe,
                                   patched_getsize,patched_NTF,*args):
    def patched_Pipe_side_effect():
      writing_side_conn,reading_side_conn = Pipe()
      patched_Pipe.wrsconn = Mock(wraps=writing_side_conn)
      patched_Pipe.rdsconn = Mock(wraps=reading_side_conn)
      return patched_Pipe.wrsconn,patched_Pipe.rdsconn
    patched_Pipe.side_effect = patched_Pipe_side_effect
    
    def patched_NTF_side_effect(*args,**kwargs):
      return_val = Mock()
      return_val.name = kwargs['prefix']
      return return_val
    patched_NTF.side_effect = patched_NTF_side_effect
    
    class TesterProc(te.multiprocessing.Process):
      def __init__(self,fifo_obj,ask_for_file_EV,shutdown):
        te.multiprocessing.Process.__init__(self)
        self.fifo_obj = fifo_obj
        self.ask_for_file = ask_for_file_EV
        self.shutdown = shutdown
      
      def run(self):
        self.fifo_obj.start_IN_end()
        while not self.ask_for_file.wait(0.01):
          pass
        self.fifo_obj.wh
        while not self.shutdown.wait(0.01):
          pass
        self.fifo_obj.close()
    
    self.fifo_obj = te.SharedFIFOfile(max_file_size_GB=0.000001,
                                      size_check_delay=0,interval_len=0.01)
    self.shutdown_EV = te.multiprocessing.Event()
    ask_for_file_EV = te.multiprocessing.Event()
    self.other_proc = TesterProc(self.fifo_obj,ask_for_file_EV,self.shutdown_EV)
    self.other_proc.start()
    self.fifo_obj.start_OUT_end()
    patched_NTF.reset_mock()
    patched_Pipe.rdsconn.reset_mock()
    
    self.assertFalse(patched_Pipe.rdsconn.recv.called)
    self.assertFalse(patched_Pipe.rdsconn.send.called)
    self.assertFalse(patched_NTF.called)
    ask_for_file_EV.set()
    
    self.shutdown_EV.set()
    self.fifo_obj.close()
    self.assertTrue(patched_NTF.called)
    patched_Pipe.rdsconn.recv.assert_called_once_with()
    patched_Pipe.rdsconn.send.assert_called_once_with('FIFOfile002_')
    self.other_proc.join()
  
  @patch('os.path.getsize')
  @patch('multiprocessing.Condition')
  def test_FIFO_close_writing_end(self,patched_Condition,patched_getsize,
                                  patched_NTF,patched_TmpDir,*args):
    def patched_Condition_side_effect(*args,**kwargs):
      patched_Condition.instance = Mock(wraps=Condition())
      return patched_Condition.instance
    patched_Condition.side_effect = patched_Condition_side_effect
    
    class TesterProc(te.multiprocessing.Process):
      def __init__(self,fifo_obj,shutdown):
        te.multiprocessing.Process.__init__(self)
        self.fifo_obj = fifo_obj
        self.shutdown = shutdown
      
      def run(self):
        self.fifo_obj.start_OUT_end()
        while not self.shutdown.wait(0.01):
          pass
        self.fifo_obj.close()
    
    self.fifo_obj = te.SharedFIFOfile(interval_len=0.01)
    self.shutdown_EV = te.multiprocessing.Event()
    self.other_proc = TesterProc(self.fifo_obj,self.shutdown_EV)
    self.other_proc.start()
    self.assertFalse(patched_Condition.instance.acquire.called)
    
    with patch('__builtin__.open',mock_open(),create=True) as patched_open:
      self.fifo_obj.start_IN_end()
    patched_Condition.instance.acquire.assert_called_once_with()
    self.shutdown_EV.set()
    self.assertFalse(patched_open.return_value.close.called)
    self.assertFalse(patched_Condition.instance.wait.called)
    self.assertFalse(patched_TmpDir.return_value.__exit__.called)
    self.assertFalse(te.SharedFIFOfile.TMPFILE.writing_side_conn.closed)
    self.assertFalse(te.SharedFIFOfile.TMPFILE.reading_side_conn.closed)
    
    self.fifo_obj.close()
    patched_open.return_value.close.assert_called_once_with()
    patched_Condition.instance.wait.assert_called_once_with(30)
    patched_TmpDir.return_value.__exit__.assert_called_once_with(None,
                                                                 None,None)
    self.assertTrue(te.SharedFIFOfile.TMPFILE.writing_side_conn.closed)
    self.assertTrue(te.SharedFIFOfile.TMPFILE.reading_side_conn.closed)
    self.other_proc.join()
  
  @patch('os.path.getsize')
  @patch('multiprocessing.Condition')
  def test_FIFO_close_reading_end(self,patched_Condition,patched_getsize,
                                  patched_NTF,*args):
    def patched_Condition_side_effect(*args,**kwargs):
      patched_Condition.instance = Mock(wraps=Condition())
      return patched_Condition.instance
    patched_Condition.side_effect = patched_Condition_side_effect
    
    def patched_NTF_ret_val_close_side_effect():
      self.assertFalse(patched_Condition.instance.notify.called)
      self.assertFalse(patched_Condition.instance.release.called)
    patched_NTF.return_value.close.side_effect = \
                                          patched_NTF_ret_val_close_side_effect
    
    class TesterProc(te.multiprocessing.Process):
      def __init__(self,fifo_obj,shutdown):
        te.multiprocessing.Process.__init__(self)
        self.fifo_obj = fifo_obj
        self.shutdown = shutdown
      
      def run(self):
        with patch('__builtin__.open',mock_open(),create=True) as patched_open:
          self.fifo_obj.start_IN_end()
        while not self.shutdown.wait(0.01):
          pass
        self.fifo_obj.close()
    
    self.fifo_obj = te.SharedFIFOfile(interval_len=0.01)
    self.shutdown_EV = te.multiprocessing.Event()
    self.other_proc = TesterProc(self.fifo_obj,self.shutdown_EV)
    self.other_proc.start()
    self.shutdown_EV.set()
    self.fifo_obj.start_OUT_end()
    self.assertTrue(self.fifo_obj.spooler.is_alive())
    self.assertFalse(patched_Condition.instance.acquire.called)
    self.assertFalse(patched_NTF.return_value.close.called)
    
    self.fifo_obj.close()
    self.other_proc.join()
    self.assertFalse(self.fifo_obj.spooler.is_alive())
    patched_Condition.instance.acquire.assert_called_once_with()
    patched_NTF.return_value.close.assert_called_once_with()
    patched_Condition.instance.notify.assert_called_once_with()
    patched_Condition.instance.release.assert_called_once_with()
  
  def tearDown(self):
    if hasattr(self,'fifo_obj'):
      if hasattr(self.fifo_obj,'spooler'):
        if self.fifo_obj.spooler.is_alive():
          self.fifo_obj.spooler.stop.set()
          self.fifo_obj.spooler.join()
    
    if hasattr(self,'shutdown_EV'):
      if not self.shutdown_EV.is_set():
        self.shutdown_EV.set()
    
    if hasattr(self,'other_proc'):
      if self.other_proc.is_alive():
        self.other_proc.join(timeout=0.1)
        if self.other_proc.is_alive():
          self.other_proc.terminate()
    
    if hasattr(te.SharedFIFOfile.TMPFILE,'writing_side_conn'):
      te.SharedFIFOfile.TMPFILE.writing_side_conn.close()
    if hasattr(te.SharedFIFOfile.TMPFILE,'reading_side_conn'):
      te.SharedFIFOfile.TMPFILE.reading_side_conn.close()


class TestSharedFIFOfileIntegration(unittest.TestCase):
  
  def setUp(self):
    reload(te)
    self.fifo_obj = te.SharedFIFOfile(top_path=None,
                                      max_file_size_GB=1.7695128917694092e-08,
                                      size_check_delay=0,interval_len=0.01)
    self.shutdown_EV = te.multiprocessing.Event()
  
  def test_writing_side_startup_rollover_and_closing(self):
    class TesterProc(te.multiprocessing.Process):
      def __init__(self,fifo_obj,call_pop,conn,shutdown):
        te.multiprocessing.Process.__init__(self)
        self.fifo_obj = fifo_obj
        self.call_pop = call_pop
        self.conn = conn
        self.shutdown = shutdown
        
      def run(self):
        self.fifo_obj.start_OUT_end()
        while not self.shutdown.is_set():
          if self.call_pop.wait(0.01):
            self.conn.send(self.fifo_obj.pop())
            self.call_pop.clear()
          else:
            continue
        self.fifo_obj.close()
      
    call_pop_EV = te.multiprocessing.Event()
    conn1,conn2 = te.multiprocessing.Pipe()
    self.other_proc = TesterProc(self.fifo_obj,call_pop_EV,conn2,
                                 self.shutdown_EV)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),0)
      
    self.other_proc.start()
    time.sleep(0.01) # Other proc needs time to run its startup
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.fifo_obj.start_IN_end()
    self.assertEqual(os.path.basename(self.fifo_obj.current_writing_file.name),
                     os.listdir(self.fifo_obj.tmpdir_obj.name)[0])
    self.fifo_obj.push(1)
    self.fifo_obj.push(2)
    self.fifo_obj.push(3)
    self.fifo_obj.push(4)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.fifo_obj.push(5)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),2)
    self.assertNotEqual(os.path.basename(self.fifo_obj.current_writing_file.name),
                        sorted(os.listdir(self.fifo_obj.tmpdir_obj.name))[0])
    self.fifo_obj.push(6)
    self.fifo_obj.push(7)
    self.fifo_obj.push(8)
    self.fifo_obj.push(9)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),3)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),1)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),2)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),3)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),4)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),3)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),5)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),2)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),6)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),7)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),8)
    call_pop_EV.set()
    self.assertEqual(conn1.recv(),9)
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.assertEqual(os.path.basename(self.fifo_obj.current_writing_file.name),
                     os.listdir(self.fifo_obj.tmpdir_obj.name)[0])
     
    self.shutdown_EV.set()
    self.fifo_obj.close()
    self.assertFalse(os.path.exists(self.fifo_obj.tmpdir_obj.name))
    self.other_proc.join()
  
  def test_reading_side_startup_rollover_and_closing(self):
    class TesterProc(te.multiprocessing.Process):
      def __init__(self,fifo_obj,call_push,conn,shutdown):
        te.multiprocessing.Process.__init__(self)
        self.fifo_obj = fifo_obj
        self.call_push = call_push
        self.conn = conn
        self.shutdown = shutdown
         
      def run(self):
        self.fifo_obj.start_IN_end()
        counter = 0
        while not self.shutdown.is_set():
          if self.call_push.wait(0.01):
            counter += 1
            self.fifo_obj.push(counter)
            self.conn.send(counter)
            self.call_push.clear()
          else:
            continue
        self.fifo_obj.close()
     
    call_push_EV = te.multiprocessing.Event()
    conn1,conn2 = te.multiprocessing.Pipe()
    self.other_proc = TesterProc(self.fifo_obj,call_push_EV,conn2,
                                 self.shutdown_EV)
    self.other_proc.start()
    self.fifo_obj.start_OUT_end()
    self.assertEqual(len(os.listdir(self.fifo_obj.tmpdir_obj.name)),1)
    self.assertEqual(os.path.basename(self.fifo_obj.current_reading_file.name),
                     os.listdir(self.fifo_obj.tmpdir_obj.name)[0])
       
    self.shutdown_EV.set()
    self.fifo_obj.close()
    self.other_proc.join()
    self.assertFalse(os.path.exists(self.fifo_obj.tmpdir_obj.name))
    # Currently the SharedFIFOfile class expects the process holding the writing
    # end to finish AFTER the process holding the reading end. When the reading
    # end's copy of tmpdir_obj is destroyed, its _closed flag is set to True
    # because the writing end is responsible for calling tmpdir_obj.__exit__().
    # The writing end's copy of tmpdir_obj is destroyed after the writing end
    # calls __exit__(). However, for the purposes of this test we want the main
    # process to hold the reading end, meaning the reading end's copy of
    # tmpdir_obj is destroyed after the writing end has called __exit__() on
    # its own copy, deleting the temp directory. However, the reading end's copy
    # does not know this, and it's _closed flag is still set to False. When the
    # reading end's copy is the destroyed on exit, tmpdir_obj.__del__() calls
    # its cleanup method and tries to delete a non-existent directory, resulting
    # in an OSError. To avoid this, we set the _closed flag to True by hand.
    self.fifo_obj.tmpdir_obj._closed = True
  
  def tearDown(self):
    if not self.shutdown_EV.is_set():
      self.shutdown_EV.set()
    
    if hasattr(self,'other_proc'):
      if self.other_proc.is_alive():
        self.other_proc.join(timeout=0.1)
        if self.other_proc.is_alive():
          self.other_proc.terminate()
    
    if os.path.exists(self.fifo_obj.tmpdir_obj.name):
      self.fifo_obj.tmpdir_obj.__exit__(None,None,None)


if __name__ == "__main__":
  #import sys;sys.argv = ['', 'Test.testName']
  unittest.main()
