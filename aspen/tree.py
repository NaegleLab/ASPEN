# Author: Roman Sloutsky <sloutsky@wustl.edu>

import copy
import weakref
from collections import Counter
from functools import wraps
from StringIO import StringIO
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

class T_BASE(object):
  
  def __new__(cls,*args,**kwargs):
    '''
    Catch cases when object passed to be wrapped is an instance of the wrapper and return
    that object without creating a new instance.
    
    Unfortunately __init__ will still be called on the existing instance, since the caller
    to __new__ (presumably type.__call__) does not see a difference between the already
    initialized object and a new instance, returned by object.__new__. To deal with this
    __init__ has to catch when this happens and not proceed with initialization, which would
    fail.
    '''
    if 'source' in kwargs:
      to_be_wrapped = kwargs['source']
    elif args:
      to_be_wrapped = args[0]
    else:
      to_be_wrapped = None
     
    if to_be_wrapped is not None and hasattr(to_be_wrapped,'wrapped'):
      return to_be_wrapped
    
    return object.__new__(cls,*args,**kwargs)
  
  def split_paths(self,paths):
    while len({p[0] for p in paths}) == 1:
      for p in paths:
        p.pop(0)
    unique_pathheads = {p[0] for p in paths}
    assert len(unique_pathheads) == 2
    new_paths = [filter(lambda x: x[0] is unique_p,paths) for unique_p in unique_pathheads]
    return [[[pset[0].pop()]] if len(pset) == 1 else pset for pset in new_paths]
  
  def recursively_split(self,paths,first_call=False,default_brlen=0.01):
    if hasattr(paths,'root'):
      return paths
        
    elif all(hasattr(c,'root') for c in paths) and len(paths) == 2:
      return Phylo.BaseTree.Clade(branch_length=default_brlen,
                                  clades=[copy.deepcopy(c.wrapped) if hasattr(c,'wrapped')
                                          else copy.deepcopy(c) for c in paths])
        
    elif all(hasattr(p,'root') or (all(hasattr(c,'root') for c in p) and len(p) > 1)
             for p in paths) and all(p[-1].name and ((hasattr(p[-2],'leaf_names')
                                                      and p[-1].name in p[-2].leaf_names)
                                                     or (p[-1].name in [l.name
                                                                        for l
                                                                        in p[-2].get_terminals()]))
                                     for p in paths if not hasattr(p,'root')):
      return self.split_paths(paths)
    else:
      if first_call:
        return self.split_paths(paths)
      else:
        return [pset.wrapped if hasattr(pset,'wrapped') else pset if hasattr(pset,'root')
                else self.recursively_split(pset) if len(pset) > 1 else pset.pop() for pset in paths]
  
  def _fix_symbols_in_names(self):
    if any([True if l.name.find(':') > -1 or l.name.find('.') > -1 else False
            for l in self._wrapped_obj.get_terminals()]):
      self._orig_leaf_names = [l.name for l in self._wrapped_obj.get_terminals()]
      for l in self._wrapped_obj.get_terminals():
        l.name = l.name.replace(':','_').replace('.','_')
  
  def _spawn(self,what_to_spawn,host=None,format=None):
    if not host:
      host = self._hostObj
    if not format:
      format = self._format
    
    spawned = T(what_to_spawn,host,format)
    spawned._spawned = {}
    spawned._spawned[id(self._wrapped_obj)] = self
    return spawned
  
  def _create_callable_attribute_access_wrapper(self,attribute):
    @wraps(getattr(self.wrapped,attribute))
    def callable_attribute_access_wrapper(*args,**kwargs):
      if attribute in ['common_ancestor','trace']:
        passed_args = [n for n in args] if len(args) > 1 else [n for n in args[0]]
        args = [n.wrapped if hasattr(n,'wrapped') else n for n in passed_args]
      elif attribute == 'is_monophyletic':
        passed_args = [n for n in args] if len(args) > 1 else [n for n in args[0]]
        args = [self.node(n).wrapped if isinstance(n,str) else n.wrapped
                if hasattr(n,'wrapped') else n for n in passed_args]
      elif attribute == 'root_with_outgroup':
        args = tuple([args[0].wrapped]+list(args[2:])) if hasattr(args[0],'wrapped') else args
      
      result = getattr(self.wrapped,attribute).__call__(*args,**kwargs)
      
      if attribute in ['get_terminals','get_nonterminals','find_clades','trace']:
        return [self._get_spawned(node) for node in result]
      elif attribute in ['common_ancestor','is_monophyletic']:
        return self._get_spawned(result) if result else result
      else:
        return result
    return callable_attribute_access_wrapper
      
  
  def __init__(self,source,host=None,format='phyloxml',keep_names_as_are=False):
    if hasattr(self,'_wrapped_obj'):
      return
    
    if isinstance(source,str) or hasattr(source,'closed'):
      try:
        self._wrapped_obj = Phylo.read(source,format)
      except:
        if format == 'phyloxml':
          if not hasattr(source,'closed') or not source.closed:
            if hasattr(source,'closed'):
              source.seek(0)
            format = 'newick'
            self._wrapped_obj = Phylo.read(source,format)
        else:
          raise
    elif hasattr(source,'count_terminals'):
      self._wrapped_obj = source
    
    self._format = format
    if not keep_names_as_are:
      self._fix_symbols_in_names()
    
    if hasattr(self._wrapped_obj,'randomized'):
      self._base_type = 'Bio.Phylo.BaseTree.Tree'
      if hasattr(self._wrapped_obj,'to_alignment'):
        self._type = 'Bio.Phylo.PhyloXML.Phylogeny'
      else:
        self._type = self._base_type
    elif hasattr(self._wrapped_obj,'is_terminal'):
      self._base_type = 'Bio.Phylo.BaseTree.Clade'
      if hasattr(self._wrapped_obj,'taxonomy'):
        self._type = 'Bio.Phylo.PhyloXML.Clade'
      else:
        self._type = self._base_type
    
    if host:
      self._hostObj = host
    elif host is None and self._base_type == 'Bio.Phylo.BaseTree.Tree':
      self._hostObj = self._wrapped_obj
    else:
      self._hostObj = 'orphan'
  
  def _get_spawned(self,item,host=None,format=None):
    if not hasattr(self,'_spawned') or self._spawned is None:
      self._spawned = {}
    if id(item) not in self._spawned:
      self._spawned[id(item)] = self._spawn(item,host,format)
    return self._spawned[id(item)]
  
  @property
  def wrapped(self):
    return self._wrapped_obj
  
  @property
  def format(self):
    return self._format
  
  @property
  def host(self):
    if self._hostObj == 'orphan':
      return self._hostObj
    elif self._hostObj is self._wrapped_obj:
      return self
    else:
      return self._get_spawned(self._hostObj)
  
  @property
  def leaf_names(self):
    return self.leaf_attrs(attr_name='name')
  
  def __deepcopy__(self,memodict={}):
    return T(copy.deepcopy(self._wrapped_obj),format=self._format,
             keep_names_as_are=True)
  
  def __getattr__(self,attribute):
    if hasattr(self,'_wrapped_obj'):
      if hasattr(self.wrapped,attribute) and hasattr(getattr(self.wrapped,attribute),'__call__'):
        return self._create_callable_attribute_access_wrapper(attribute)
      else:
        if attribute == 'root':
          return self._get_spawned(self.wrapped.root)
        elif attribute == 'clades':
          call_on = self if hasattr(self.wrapped,'clades') else self.root
          return [call_on._get_spawned(clade) for clade in call_on.wrapped.clades]
        elif not hasattr(self.wrapped,attribute) and\
                                          self._base_type == 'Bio.Phylo.BaseTree.Tree':
          return getattr(self.wrapped.root,attribute)
        else:
          return getattr(self.wrapped,attribute)
    else:
      raise AttributeError("T instance has no attribute _wrapped_obj")
  
  def __setattr__(self,attribute,value):
    if attribute == 'clades':
      call_on = self.wrapped if hasattr(self.wrapped,'clades') else self.wrapped.root
      call_on.clades = [n.wrapped if hasattr(n,'wrapped') else n for n in value]
    elif hasattr(self,'_wrapped_obj') and hasattr(self.wrapped,attribute):
      setattr(self._wrapped_obj,attribute,value)
    elif hasattr(self,'_wrapped_obj') and hasattr(self.wrapped.root,attribute):
      setattr(self._wrapped_obj.root,attribute,value)
    else:
      super(T_BASE,self).__setattr__(attribute,value)
  
  def __getstate__(self):
    return {'wrapped':self._wrapped_obj,'format':self._format}
  
  def __setstate__(self,state_dict):
    self.__init__(source=state_dict['wrapped'],format=state_dict['format'],
                  keep_names_as_are=True)
  
  def __eq__(self,other):
    if hasattr(self,'nested_set_repr') and hasattr(other,'nested_set_repr'):
      return self.nested_set_repr() == other.nested_set_repr()
    else:
      return False
  
  def __ne__(self,other):
    return not self.__eq__(other)
  
  def __repr__(self):
    return self.__class__.__name__+'( '+repr(self.wrapped)+' )'
  
  def leaf_attrs(self,attr_name='name'):
    def handle_returning(return_this):
      if hasattr(return_this,'pop') and len(return_this) == 1:
        return return_this.pop()
      else:
        return return_this
    if attr_name == 'name':
      return handle_returning([l.name for l in self.wrapped.get_terminals()])
    elif attr_name == 'sciname':
      return handle_returning([l.taxonomy.scientific_name for l in self.wrapped.get_terminals()])
    elif attr_name == 'species':
      return handle_returning([l.properties[0].value for l in self.wrapped.get_terminals()])
  
  def node(self,how_to_find,leaf_count=None):
    if isinstance(how_to_find,float):
      matches = [self._get_spawned(n) for n in self.wrapped.get_nonterminals()
                 if n.branch_length == how_to_find]
      if len(matches) == 0:
        print "No internal nodes with branch length",how_to_find
        return None
      if len(matches) > 1:
        print "Multiple internal nodes with branch length",how_to_find
        print "Cannot perform leaf count check"
        print "Returning all matches"
        return matches
      elif leaf_count is not None:
        assert [n for n in self.wrapped.get_nonterminals() if
                n.branch_length == how_to_find].pop().count_terminals() == leaf_count
      return matches.pop()
    elif isinstance(how_to_find,str):
      matches = [self._get_spawned(n) for n
                 in self.wrapped.find_clades(how_to_find)]
      if len(matches) == 0:
        print "No internal nodes or leaves with name", how_to_find
        return None
      if len(matches) > 1:
        print "Multiple internal nodes or leaves with name",how_to_find
        print "Returning all matches"
        return matches
      else:
        return matches.pop()
  
  def get_path(self,target=None,full=False,**kwargs):
    if full:
      path_root = self._hostObj
    else:
      path_root = self.wrapped
    if not target and not kwargs:
      return [self._get_spawned(n) for n in self._hostObj.get_path(self.wrapped)]
    elif hasattr(target,'wrapped'):
      return [self._get_spawned(n) for n in path_root.get_path(target.wrapped,**kwargs)]
    else:
      return [self._get_spawned(n) for n in path_root.get_path(target,**kwargs)]
  
  def trace_dist(self,node1,node2=None):
    if node2 is None:
      return len(self.get_path(node1))
    else:
      return len(self.trace(node1,node2))-1
  
  def nested_set_repr(self):
    # It would be great to store the computed result for fast retrieval next time it's
    # needed, but there's a catch: I've introduced mechanisms to modify trees, so the
    # representation might be out of date. Same goes for leaf names and any number of
    # computed properties that could be stored.
    # If there was a flag that the object could check to see whether it has been modified
    # since a property was last computed, it would solve this problem. I just have to make
    # sure the flag is checked in all the appropriate places. Sounds like a job for a
    # decorator ...
    return frozenset(c.name if c.is_terminal() else c.nested_set_repr() for
                     c in self.clades)
  
  @property
  def parent(self):
    if len(self.get_path()) == 1:
      if hasattr(self._hostObj,'root'):
        return self._hostObj.root
      else:
        return self._hostObj
    else:
      return self.get_path()[-2]
  
  def write(self,filepath,format=None,**kwargs):
    if not format:
      format = self.format
    if filepath == 'as_string':
      strhandle = StringIO()
      Phylo.write(self.wrapped,strhandle,format,**kwargs)
      return strhandle.getvalue()
    else:
      Phylo.write(self.wrapped,filepath,format,**kwargs)



class T(T_BASE):
  _to_wrappers_map = weakref.WeakValueDictionary()
  _to_wrapped_map = weakref.WeakValueDictionary()
  
  _keep_alive = {}
  _pickle_counts = Counter()
  
  @classmethod
  def _nsrepr(cls,obj):
    if obj.is_terminal():
      return obj.name
    elif all(c.is_terminal() for c in obj.clades):
      return frozenset(c.name for c in obj.clades)
    elif any(c.is_terminal() for c in obj.clades):
      return frozenset(c.name if c.is_terminal() else
                       cls._nsrepr(c) for c in obj.clades)
    else:
      return frozenset(cls._nsrepr(c) for c in obj.clades)
  
  @classmethod
  def _check_clade(cls,clade_obj):
    try:
      assert cls._to_wrapped_map[cls._nsrepr(clade_obj)] is clade_obj
    except KeyError:
      assert False,"Unknown non-leaf clade "+repr(cls._nsrepr(clade_obj))
    except AssertionError:
      assert False,"Duplicate clade "+repr(cls._nsrepr(clade_obj))
  
  def __new__(cls,*args,**kwargs):
    if 'source' in kwargs:
      obj = kwargs['source']
    elif args:
      obj = args[0]
    else:
      return T_BASE.__new__(cls,*args,**kwargs)
    if hasattr(obj,'is_terminal'):
      try:
        cls._check_clade(obj)
      except AssertionError:
        print "Caught clade error in __new__"
        raise
      obj_to_wrapper_key = id(obj)
      if obj_to_wrapper_key in cls._to_wrappers_map:
        return cls._to_wrappers_map[obj_to_wrapper_key]
      else:
        new_obj = T_BASE.__new__(cls,*args,**kwargs)
        cls._to_wrappers_map[obj_to_wrapper_key] = new_obj
        return new_obj
    else:
      return T_BASE.__new__(cls,*args,**kwargs)
  
  def _get_spawned(self,item,host=None,format=None):
    if not hasattr(self,'_spawned') or self._spawned is None:
      self._spawned = {}
    id_item = id(item)
    if id_item not in self._spawned:
      if not host:
        host = self._hostObj
      if not format:
        format = self._format
      self._spawned[id_item] = type(self)(item,host,format)
    return self._spawned[id_item]
  
  @classmethod
  def requisition(cls,clade_def1,clade_def2,*args):
    new_clades_attr = []
    for c in (clade_def1,clade_def2)+args:
      if isinstance(c,str):
        if c in cls._to_wrapped_map:
          new_clades_attr.append(cls._to_wrapped_map[c])
        else:
          new_leaf = Clade(name=c)
          cls._to_wrapped_map[c] = new_leaf
          new_clades_attr.append(new_leaf)
      else:
        try:
          cls._check_clade(c)
        except AssertionError:
          print "Caught clade error in requisition()"
          raise
        new_clades_attr.append(c)
    proposed_key = frozenset(cls._nsrepr(c) for c in new_clades_attr)
    if proposed_key in cls._to_wrapped_map:
      return cls(cls._to_wrapped_map[proposed_key])
    else:
      nonleaf_clade = Clade(clades=new_clades_attr)
      cls._to_wrapped_map[cls._nsrepr(nonleaf_clade)] = nonleaf_clade
      return cls(nonleaf_clade)
  
  @classmethod
  def rebuild_on_unpickle(cls,clade_repr,top_level_call=True):
    if clade_repr in cls._to_wrapped_map:
      if top_level_call:
        return cls(cls._to_wrapped_map[clade_repr])
      else:
        return cls._to_wrapped_map[clade_repr]
    elif isinstance(clade_repr,str):
      leaf = Clade(name=clade_repr)
      cls._to_wrapped_map[clade_repr] = leaf
      return leaf
    if top_level_call:
      return cls.requisition(*[cls.rebuild_on_unpickle(c_repr,False)
                               for c_repr in clade_repr])
    else:
      return cls.requisition(*[cls.rebuild_on_unpickle(c_repr,False)
                               for c_repr in clade_repr]).wrapped
  
  def check_in_pickle(self):
    key = self._nsrepr(self._wrapped_obj)
    if key not in self._keep_alive:
      assert self._pickle_counts[key] == 0
      self._keep_alive[key] = self._wrapped_obj
    self._pickle_counts[key] += 1
    return key
  
  @classmethod
  def check_out_pickle(cls,key):
    assert key in cls._pickle_counts and cls._pickle_counts[key] > 0
    assert key in cls._keep_alive
    cls._pickle_counts[key] -= 1
    if cls._pickle_counts[key] == 0:
      kept_alive = cls._keep_alive.pop(key)
      cls._pickle_counts.pop(key)
    else:
      kept_alive = cls._to_wrapped_map[key]
    return cls(kept_alive)
