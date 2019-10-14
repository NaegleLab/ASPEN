# Author: Roman Sloutsky <sloutsky@wustl.edu>

from topolenum import *
import subprocess
from sys import float_info,stderr,argv

#===============================================================================
# IO utils
#===============================================================================

def load_pwhist(fpath):
  with open(fpath) as fh:
    return [(eval(l.split('\t')[0]),eval(l.split('\t')[1])) for l in fh]

def write_enumeration_results(results,targetpath):
  with open(targetpath,'w') as wh:
    for r in results:
      wh.write(repr(r))


#===============================================================================
# Running and monitoring topology assembly by enumeration
#===============================================================================
 
class EnumerationObserver(object):
  MINUTE = 60
  HOUR = 3600
  DAY = 86400
   
  def __init__(self,terminate_after=None,terminator_file='stop_enumeration',
               timestamp_frequency=HOUR,report_frequency=None,
               username_for_top=None,report_on_workers=False):
    if isinstance(terminate_after,basestring):
      if terminate_after[-1] == 'h':
        self.terminate_after = float(terminate_after[:-1])*self.HOUR
      elif terminate_after[-1] == 'm':
        self.terminate_after = float(terminate_after[:-1])*self.MINUTE
      else:
        self.terminate_after = float(terminate_after)
    else:
      self.terminate_after = terminate_after
    self.terminator = terminator_file
    self.timestamp_freq = timestamp_frequency
    self.username = username_for_top
    self.report_freq = report_frequency
    self.report_on_workers = report_on_workers
    self.old_min_score = -float_info.max
    self.start_time = time.time()
    self.time_of_last_stamp = self.start_time
    self.time_of_last_report = self.start_time
     
    self.first_call = True
    self.interrupt_reported = False
   
  def time_since(self,time_in_past):
    return time.time()-time_in_past
   
  @property
  def total_elapsed_time(self):
    return self.time_since(self.start_time)
   
  @property
  def timestamp(self):
    stamp = ''
    seconds = self.total_elapsed_time
    true_seconds = seconds
    if true_seconds >= self.DAY:
      stamp += '%0.0fd:' %  ((seconds-seconds % self.DAY) / self.DAY)
      seconds = seconds % self.DAY
    if true_seconds >= self.HOUR:
      stamp += '%0.0fh:' %  ((seconds-seconds % self.HOUR) / self.HOUR)
      seconds = seconds % self.HOUR
    if true_seconds >= self.MINUTE:
      stamp += '%0.0fm:' %  ((seconds-seconds % self.MINUTE)/self.MINUTE)
      seconds = seconds % self.MINUTE
    stamp += '%0.0fs\t  ' % seconds
    return stamp
   
  def report_timestamp(self):
    print >>stderr,self.timestamp,"Elapsed since start"
    self.time_of_last_stamp = time.time()
  
  def report_top_output(self,enum_proc,workers=None):
    self.dict_proc_PID = enum_proc.encountered_assemblies_manager._process.pid
    self.queue_proc_PID = enum_proc.assembly_queue_manager._process.pid
    print >>stderr,'='*80
    subprocess.call('top -n 1 -b | grep PID',shell=True)
    print >>stderr,'-'*31+'Shared Dict Process'+'-'*30
    grep_this = '"'+str(self.dict_proc_PID)+' %s"' % self.username
    subprocess.call('top -n 1 -u %s -b | grep ' % self.username + grep_this,
                    shell=True)
    if workers is not None:
      print >>stderr,'-'*32+'Worker Processes'+'-'*32
      proc_PIDs = workers.keys() + [self.queue_proc_PID,self.dict_proc_PID]
      proc_args = ' '.join(['-p'+str(pid) for pid in proc_PIDs])
      subprocess.call('top -n 1 '+proc_args+' -b | grep %s' % self.username,
                      shell=True)
    print >>stderr,'~'*80
    self.time_of_last_report = time.time()
   
 
  def report_score(self,enum_proc):
    if enum_proc.min_score.value > self.old_min_score:
      self.old_min_score = enum_proc.min_score.value
      print >>stderr,self.timestamp,"New worst score is",self.old_min_score
   
  def report_interrupt(self,enum_proc):
    print >>stderr,self.timestamp,"Interrupt requested, writing save to",\
                   enum_proc.save_file_name
   
  def proceed_permission_check(self):
    if self.terminator in os.listdir('.'):
      os.remove(self.terminator)
      return False
    elif self.terminate_after is not None and\
                                self.total_elapsed_time > self.terminate_after:
      return False
    else:
      return True
     
   
  def __call__(self,enum_proc,workers,interrupt=False):
    if self.time_since(self.time_of_last_stamp) > self.timestamp_freq:
      self.report_timestamp()
    if self.report_freq is not None:
      if self.first_call or\
                  self.time_since(self.time_of_last_report) > self.report_freq:
        if self.username is not None:
          if self.report_on_workers:
            self.report_top_output(enum_proc,workers)
          else:
            self.report_top_output(enum_proc)
          self.first_call = False
    if interrupt:
      if not self.interrupt_reported:
        self.report_interrupt(enum_proc)
        self.interrupt_reported = True
    else:
      self.report_score(enum_proc)
 
 
def enumerate_topologies(leafdist_histograms,restart_from=None,
                         proceed_permission_callable=lambda: True,
                         wait_duration=10,observer_callable=None,**kwargs):
  enumeration_proc = MainTopologyEnumerationProcess(leafdist_histograms,
                                                    restart_from=restart_from,
                                                    **kwargs)
  enumeration_proc.start()
  workers = {}
  while len(workers) < enumeration_proc.num_workers:
    if enumeration_proc.get_PIDs.poll(1):
      new_worker = enumeration_proc.get_PIDs.recv()
      workers[new_worker[0]] = new_worker[1]
   
  while not enumeration_proc.finished.wait(wait_duration):
    if observer_callable is not None:
      observer_callable(enumeration_proc,workers)
    if not proceed_permission_callable():
      enumeration_proc.stop.set()
      break
   
  if enumeration_proc.stop.is_set():
    results = None
    while not enumeration_proc.save_written.wait(wait_duration):
      if observer_callable is not None:
        observer_callable(enumeration_proc,workers,interrupt=True)
  elif enumeration_proc.finished.is_set():
    results = []
    finished_worker_counter = 0
    while finished_worker_counter < \
                      enumeration_proc.expected_number_results_queue_sentinels:
      received = enumeration_proc.results_queue.get()
      if received == 'FINISHED':
        finished_worker_counter += 1
      else:
        results.append(received)
    results.sort(key=lambda x: x.score,reverse=True)
    while len(results) > enumeration_proc.num_requested_topologies:
      results.pop()
  enumeration_proc.shutdown.set()
  enumeration_proc.join()
  enumeration_proc.clean_up()
   
  return results

if __name__ == '__main__' and len(argv) > 1:
#------------------------------------------------------------------------------
# Arguments

  pwhistfile = argv[1]
  numproc = int(argv[2]) if len(argv) > 2 else 4
  
  numtopologies = 1000
  restartfrom = None
  walltime = None
  savefilename = 'early_termination_save'
  timestampfreq = 30
  distconstraintfreq = 0.99
  absdistfreq = 0.001
  maxfifofilesize = 0.001
  workspacesize = 100
  topoffacceptanceratio = 2.0
  topoffacceptancestiffness = 1.0
  outfile = 'aspen_topologies.txt'

#------------------------------------------------------------------------------
# Run

  pwhist = load_pwhist(pwhistfile)
  observer=EnumerationObserver(terminate_after=walltime,
                               timestamp_frequency=timestampfreq)
  check_proceed_permission = observer.proceed_permission_check
  num_workers = numproc
  results = enumerate_topologies(pwhist,restart_from=restartfrom,
                          proceed_permission_callable=check_proceed_permission,
                                 observer_callable=observer,
                                 num_workers=num_workers,
                                 num_requested_topologies=numtopologies,
                                 constraint_freq_cutoff=distconstraintfreq,
                                 absolute_freq_cutoff=absdistfreq,
                                 fifo_max_file_size=maxfifofilesize,
                                 max_workspace_size=workspacesize,
                                 acceptance_ratio_param=topoffacceptanceratio,
                          acceptance_stiffness_param=topoffacceptancestiffness,
                                 save_file_name=savefilename)
  if results is not None:
    write_enumeration_results(results,outfile)

