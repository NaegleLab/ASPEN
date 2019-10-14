import os
import tempfile
import threading
import multiprocessing
import cPickle as pickle
from .tempdir import TemporaryDirectory

#===============================================================================
# On-disc FIFO (for storing mid-construction topologies)
#===============================================================================


class FIFOfile(object):
  class TMPFILE(object):
    def __new__(cls,*args,**kwargs):
      cls.instcount += 1
      return object.__new__(cls,*args,**kwargs)
    
    @classmethod
    def init_class(cls,mode,wbuf,rbuf,suffix,delete,dir,check_delay):
      cls.mode = mode
      cls.wbuf = wbuf
      cls.rbuf = rbuf
      cls.suffix = suffix
      cls.delete = delete
      cls.dir = dir
      cls.check_delay = check_delay
      cls.instcount = 0
    
    @classmethod
    def start_spooling(cls):
      cls.file_spool = []
    
    @classmethod
    def spool(cls):
      cls.file_spool.append(cls())
      return cls.file_spool[-1]
    
    @classmethod
    def pop_from_spool(cls):
      if cls.file_spool:
        return cls.file_spool.pop(0)
      else:
        return None
    
    def __init__(self):
      self.rh = tempfile.NamedTemporaryFile('r'+self.mode,self.rbuf,
                                            suffix=self.suffix,
                                            prefix='FIFOfile'+\
                                            str(self.instcount).zfill(3)+'_',
                                            dir=self.dir,delete=self.delete)
      self.name = self.rh.name
      self._size = 0.0
      self.access_count_since_size_check = 0
    
    @property
    def size(self):
      if self.access_count_since_size_check >= self.check_delay:
        self._size = os.path.getsize(self.name)
        self.access_count_since_size_check = 0
      else:
        self.access_count_since_size_check += 1
      return self._size
    
    def open(self):
      self.wh = open(self.name,'w'+self.mode,self.wbuf)
    
    def close(self):
      self.wh.close()
    
    def discard(self):
      self.rh.close()
      assert os.path.exists(self.name) == False
    
  
  def __init__(self,mode='b',wbuffering=0,rbuffering=0,delete=True,top_path='.',
               suffix='',max_file_size_GB=1.0,size_check_delay=100):
    self.tmpdir_obj = TemporaryDirectory(dir=top_path,prefix='FIFOworkspace_',
                                         suffix=suffix)
    self.max_size = max_file_size_GB*1024**3 #os.path.getsize() reports in bytes
    self.TMPFILE.init_class(mode,wbuffering,rbuffering,suffix,delete,
                            os.path.realpath(self.tmpdir_obj.__enter__()),
                            size_check_delay)
    
  def start_OUT_end(self):
    self.current_reading_file = self.TMPFILE()
    self.TMPFILE.start_spooling()
  
  def start_IN_end(self):
    self.current_writing_file = self.current_reading_file
    self.current_writing_file.open()
  
  @property
  def wh(self):
    if self.current_writing_file.size > self.max_size:
      self.current_writing_file.close()
      self.current_writing_file = self.TMPFILE.spool()
      self.current_writing_file.open()
    return self.current_writing_file.wh
  
  @property
  def rh(self):
    old_pos = self.current_reading_file.rh.tell()
    self.current_reading_file.rh.read(1)
    if self.current_reading_file.rh.tell() == old_pos: # => cursor is at EOF
      next_file = self.TMPFILE.pop_from_spool()
      # We can't check whether the writing handle is closed because it may not
      # exist in this process. However, if another file is already spooled,
      # this file must be closed and safe to discard.
      if next_file is not None:
        self.current_reading_file.discard()
        self.current_reading_file = next_file
        # Need to update old_pos after rolling over
        old_pos = self.current_reading_file.rh.tell()
    # If the reading operation successfully advanced the cursor, put it back!
    # If it didn't, we still want to do the same thing:
    #
    # On some platforms performing a read operation when at EOF appears to
    # put the file handle into a state where further read operations return
    # nothing and fail to advance the cursor, even after additional writing.
    # When in this state, calling handle.seek(pos), where pos is the output
    # of handle.tell() seems to revert the handle to a normal reading state,
    # fixing the problem.
    #
    # If rollover was successful, this will leave the cursor at the top
    self.current_reading_file.rh.seek(old_pos)
    return self.current_reading_file.rh
  
  def pop(self):
    try:
      result = pickle.load(self.rh)
    except EOFError:
      result = None
    return result
  
  def push(self,item):
    pickle.dump(item,self.wh,pickle.HIGHEST_PROTOCOL)
  
  def close(self):
    if self.current_reading_file is self.current_writing_file:
      assert not self.TMPFILE.file_spool
      self.current_writing_file.close()
      self.current_reading_file.discard()
    else:
      assert self.current_reading_file.wh.closed
      self.current_reading_file.discard()
      for tmpfile in self.TMPFILE.file_spool:
        tmpfile.close()
        tmpfile.discard()
    self.tmpdir_obj.__exit__(None,None,None)


class SharedFIFOfile(FIFOfile):
  class TMPFILE(FIFOfile.TMPFILE):
    @classmethod
    def init_class(cls,*args,**kwargs):
      super(SharedFIFOfile.TMPFILE,cls).init_class(*args,**kwargs)
      cls.writing_side_conn,cls.reading_side_conn = multiprocessing.Pipe()
    
    @classmethod
    def spool(cls):
      cls.writing_side_conn.send(None)
      return cls(cls.writing_side_conn.recv())
    
    def __init__(self,filename=None):
      if filename is None:
        FIFOfile.TMPFILE.__init__(self)
      else:
        self.name = filename
        self._size = os.path.getsize(filename)
        self.access_count_since_size_check = 0
  
  
  class SpoolerThread(threading.Thread):
    def __init__(self,connection,spool_callable,interval_len=1):
      proc_name = multiprocessing.current_process().name
      threading.Thread.__init__(self,name='--'.join([proc_name,
                                                     'TmpfileSpoolerThread']))
      self.stop = threading.Event()
      self.connection = connection
      self.spool = spool_callable
      self.interval_len = interval_len
    
    def run(self):
      while not self.stop.is_set():
        if self.connection.poll(self.interval_len):
          assert self.connection.recv() is None
          self.connection.send(self.spool().name)
  
  
  def __init__(self,*args,**kwargs):
    if 'interval_len' in kwargs:
      self.spooler_polling_interval_len = kwargs.pop('interval_len')
    else:
      self.spooler_polling_interval_len = 1
    FIFOfile.__init__(self,*args,**kwargs)
    
    self.lock = multiprocessing.Lock()
    self.acquire = self.lock.acquire
    self.release = self.lock.release
    
    self.event = multiprocessing.Event()
    self.set = self.event.set
    self.is_set = self.event.is_set
    self.clear = self.event.clear
    self.wait = self.event.wait
    
    self.shutdown_baton = multiprocessing.Condition()
  
  def start_OUT_end(self):
    self.side = 'reading'
    FIFOfile.start_OUT_end(self)
    self.TMPFILE.reading_side_conn.send(self.current_reading_file.name)
    self.spooler = self.SpoolerThread(self.TMPFILE.reading_side_conn,
                                      super(self.TMPFILE,self.TMPFILE).spool,
                                      self.spooler_polling_interval_len)
    self.spooler.start()
  
  def start_IN_end(self):
    self.side = 'writing'
    self.current_writing_file = self.TMPFILE(self.TMPFILE.writing_side_conn.recv())
    self.current_writing_file.open()
    self.shutdown_baton.acquire()
  
  def _sync_safe_method_call(self,method,args,already_have_lock=False):
    if not already_have_lock:
      with self.lock:
        return method(self,*args)
    else:
      return method(self,*args)
  
  def pop(self):
    self.wait()
    result = self._sync_safe_method_call(FIFOfile.pop,tuple())
    if result is None:
      self.clear()
    return result
  
  def push(self,item,already_have_lock=False):
    self._sync_safe_method_call(FIFOfile.push,(item,),already_have_lock)
    if not already_have_lock:
      self.set()
  
  def push_all(self,items):
    with self.lock:
      for item in items:
        self.push(item,already_have_lock=True)
    self.set()
  
  def close(self):
    if self.side == 'reading':
      self.spooler.stop.set()
      self.spooler.join(timeout=10)
      self.shutdown_baton.acquire()
      self.current_reading_file.discard()
      for tmpfile in self.TMPFILE.file_spool:
        tmpfile.discard()
      self.shutdown_baton.notify()
      self.shutdown_baton.release()
    else:
      self.current_writing_file.close()
      self.set() # Free QueueLoader proc from wait in pop() so it can shut down
      self.shutdown_baton.wait(30)
      self.tmpdir_obj.__exit__(None,None,None)
      self.TMPFILE.writing_side_conn.close()
      self.TMPFILE.reading_side_conn.close()
