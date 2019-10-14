# Author: Roman Sloutsky <sloutsky@wustl.edu>

import sys
import time
import os
import math
import itertools
import multiprocessing
import Queue
import gc
import tarfile
from cStringIO import StringIO
from collections import defaultdict,namedtuple,Hashable
from .tree import T_BASE,T
from fifo import FIFOfile,SharedFIFOfile

#===============================================================================
# Topology assembly extension through branching
#===============================================================================


LPDF = namedtuple('LeafPairDistanceFrequency',['leaves','dist','freq'])


class ProposedExtension(object):
  
  IndexedPair = namedtuple('IndexedPair',['index','pair'])
  
  def __init__(self,child1,child2):
    if any([isinstance(child1,str),isinstance(child2,str)]):
      # This extension is an attachment of a new leaf to a built clade ...
      # ... so there should be one of each.
      assert isinstance(child1,str) != isinstance(child2,str)
      self.new_leaf = child1 if isinstance(child1,str) else child2 # Figure out ...
      self.built_clade = child2 if isinstance(child1,str) else child1 # .. which is which
      # When a new leaf is attached to a built clade via a new root, the distance between
      # the new leaf and each of the leaves in the built clade will be distance of leaf
      # in existing clade to its current root + 1 to account for the new root
      self.unverified = dict((frozenset({leaf,self.new_leaf}),
                                     self.built_clade.clade.trace_dist(leaf)+1
                              )
                                    for leaf in self.built_clade.clade.leaf_names
                             )
    else: # This extension is the joining of two built clades
      self.clades = sorted([child1,child2],key=lambda x: x.index,reverse=True)
      # When two built clades are joined, the resulting distance between any pair of leaves
      # such that each leaf belongs to a different clade will be
      # sum(distance of each leaf to its current root) + 1 to account for the new root.
      self.unverified = dict((frozenset({leafpair[0][0],leafpair[1][0]}),
                                              leafpair[0][1]+leafpair[1][1]+1
                                     )
                                    for leafpair
                                    in itertools.product(*([(leaf,
                                                             clade.clade.trace_dist(leaf)
                                                             )
                                                            for leaf in clade.clade.leaf_names
                                                            ]
                                                           for clade in self.clades
                                                           )
                                                         )
                                    )
    self.consistent = {}
    self.inconsistent = {}
    self.verified = set()
    self.score = 0.0
  
  def check_pair(self,pair,i):
    if pair.leaves in self.consistent:
      # If this leaf pair has already been added to consistent, then ...
      # ... it should not be in unverified any more ...
      assert pair.leaves not in self.unverified
      # ... and this new distance should be different from the consistent one ...
      assert pair.dist != self.consistent[pair.leaves].pair.dist
      self.inconsistent[i] = pair # ... so it goes into inconsistent
    else: # If it hasn't been added to consistent ...
      if pair.dist == self.unverified[pair.leaves]: # ... and its distance matches expected
        self.consistent[pair.leaves] = self.IndexedPair(i,pair) # it goes into consistent
        # and is "verified", so pop it from unverified and add to verified
        self.unverified.pop(pair.leaves)
        self.verified.add(pair.leaves)
        self.score += math.log(pair.freq)
      else: # ... and its distance doesn't match expected
        self.inconsistent[i] = pair # it goes into inconsistent ...
        # ... and it remains "unverified", so don't pop it from unverified
  
  def build_extension(self,assemblyobj,in_place=False):
    if not in_place:
      assemblyobj = assemblyobj.copy()
    assert not self.unverified
    
    # Assemble for removal from constraints_idx indeces of consistent and inconsistent pairdists
    drop_these = [pair.index for pair in self.consistent.values()]+self.inconsistent.keys()
    
    # Make sure to use clade(s) from this assemblyobj for construction, not the ones from
    # the original that are stored in self.clades
    if hasattr(self,'new_leaf'):
      built_clade = assemblyobj.built_clades.pop(self.built_clade.index)
      assemblyobj.free_leaves.remove(self.new_leaf)
      assert set(self.built_clade.clade.leaf_names) == set(built_clade.leaf_names)
      new_clades_attr = [built_clade.wrapped,self.new_leaf]
      # One more thing to do if this extension is an attachment of a new leaf:
      # Add for removal constraint pairdists where the new leaf has distance 1 with any
      # other leaf, regardless of whether the other leaf is involved in this extension.
      # All such constraint distances are inconsistent with this extension. This additional
      # check is necessary because of the special way new pair extensions are handled - with
      # no questions asked.
      drop_these.extend([i for i in assemblyobj.constraints_idx
                         if self.new_leaf in assemblyobj.constraints_master[i].leaves and
                                             assemblyobj.constraints_master[i].dist == 1])
    else:
      assert all(set(c.clade.leaf_names)==set(assemblyobj.built_clades[c.index].leaf_names)
                                                            for c in self.clades)
      new_clades_attr = [assemblyobj.built_clades.pop(c.index).wrapped for c in self.clades]
    
    assemblyobj.constraints_idx = filter(lambda x: x not in drop_these,
                                         assemblyobj.constraints_idx)
    
    assemblyobj.built_clades.append(T.requisition(*new_clades_attr))
    assemblyobj.recompute(extension=self)
    assemblyobj.score += self.score
    return assemblyobj


class TreeAssembly(object):
  
  IndexedClade = namedtuple('IndexedClade',['index','clade'])
  
  class KeyPassingDefaultDict(defaultdict):
    def __missing__(self,key):
      self[key] = self.default_factory(key)
      return self[key]
  
  def __init__(self,pwleafdist_histograms,constraint_freq_cutoff,leaves_to_assemble,
               absolute_freq_cutoff=0.01,keep_alive_when_pickling=True):
    #===========================================================================
    # The data attributes below will are set on the class, not in instances,
    # meaning they will be shared between all instances of this class, saving
    # lots of space as the number of assemblies grows.
    # *** The are INVARIANTS and should not be changed by instances!!! ***
    #===========================================================================
    
    # Make a sorted master reference tuple of pw distance freq constraints, ...
    type(self).constraints_master = tuple(sorted([
                                            # by creating LPDF tuples ...
                                            LPDF(pair_histogram[0],score[0],score[1])
                                            # for every histogram of pw distances for a leaf pair ...
                                            for pair_histogram in
                                            # present in a built-on-the-fly subset of
                                            # the full histogram for that pair ...
        [(leafpair,[dist for i,dist in enumerate(distances) if
                    # of distances ordered by freq such that the sum of their freqs
                    # is less than the requested frequency cutoff, ...
                    sum(d[1] for d in distances[:i]) < constraint_freq_cutoff]
          ) for leafpair,distances in pwleafdist_histograms]
                                            for score in pair_histogram[1]],
                                           # sorted on (shortest dist, highest freq)
                                           key=lambda x: (x.dist,1-x.freq)
                                                 )
                                          )
    
    # Unlike constraints_master, which guarantees invariance by being recursively immutable,
    # this is a mutable dict. Is there a way to force immutability?
    type(self).pwdist_histograms_dict = {leafpair:dict(dist_histogram)
                                         for leafpair,dist_histogram in pwleafdist_histograms}
    
    # More technically mutable variables that should not be changed by instances
    type(self).abs_cutoff = absolute_freq_cutoff
    type(self).leaves_master = set(leaves_to_assemble)
    type(self).pickle_encoding = {chr(i):t for i,t in enumerate('[],;')}
    type(self).pickle_encoding.update({chr(i+4):repr(l) for i,l in
                                       enumerate(type(self).leaves_master)})
    type(self).total_nodes_to_build = len(leaves_to_assemble) - 1
    type(self).best_possible = sum(math.log(max(p[1],key=lambda x: x[1])[1])
                                   for p in pwleafdist_histograms)
    type(self).keep_alive = keep_alive_when_pickling
    
    #===========================================================================
    # END of class attributes
    #===========================================================================
    
    self.built_clades = []
    self.free_leaves = set(leaves_to_assemble)
    self.constraints_idx = range(len(self.constraints_master))
    self.score = 0.0
  
  def rebuild_constraints_idx(self):
    self.constraints_idx = [i for i,_ in enumerate(self.constraints_master)]
    intra_clade_pairs = [frozenset(p) for c in self.built_clades for p in
                                        itertools.combinations(c.leaf_names,2)]
    drop_these_idx = [i for i,d in enumerate(self.constraints_master)
                                        if d.leaves in intra_clade_pairs][::-1]
    for i in drop_these_idx:
      self.constraints_idx.pop(i)
    l_accounted_for = set.union(*[set(leafset) for leafset in
                                                     self.pairs_accounted_for])
    drop_these_idx = [idx for idx,i in enumerate(self.constraints_idx)
                      if self.constraints_master[i].dist == 1 and
                      self.constraints_master[i].leaves & l_accounted_for][::-1]
    for i in drop_these_idx:
      self.constraints_idx.pop(i)
  
  def convert_containers(self,convert_this,container_type=None):
    result = []
    for member in convert_this:
      if isinstance(member,Hashable) and member in self.leaves_master:
        result.append(member)
      else:
        result.append(self.convert_containers(member,container_type))
    if container_type is not None:
      return container_type(result)
    else:
      return result
  
  def __getstate__(self):
    state = {'score':self.score}
    if self.keep_alive:
      state['built_clades'] = [c.check_in_pickle() for c in self.built_clades]
    else:
      state['built_clades'] = [T._nsrepr(c) for c in self.built_clades]
    state['built_clades'] = ';'.join(''.join(repr(
                                          self.convert_containers(c)).split())
                                     for c in state['built_clades'])
    for k,v in self.pickle_encoding.items():
      state['built_clades'] = state['built_clades'].replace(v,k)
    return state['built_clades'],state['score'],self.best_case,\
                                                       self.nodes_left_to_build
  
  def _unpack_state(self,state):
    s = StringIO(state[0])
    clades = []
    while True:
      read_byte = s.read(1)
      if not read_byte:
        break
      else:
        clades.append(self.pickle_encoding[read_byte])
    return {'score':state[1],'_best_case':state[2],'_nodes_left_to_build':state[3],
            'built_clades':[self.convert_containers(eval(m),frozenset)
                            for m in ''.join(clades).split(';')]}
  
  def __setstate__(self,state):
    state = self._unpack_state(state)
    for k,v in state.items():
      if k != 'built_clades':
        self.__dict__[k] = v
    if self.keep_alive:
      self.__dict__['built_clades'] = [T.check_out_pickle(k)
                                       for k in state['built_clades']]
    else:
      self.__dict__['built_clades'] = [T.rebuild_on_unpickle(clade_repr)
                                       for clade_repr in state['built_clades']]
    self.free_leaves = self.leaves_master -\
                          set.union(*[set(leafset) for leafset in
                                                     self.pairs_accounted_for])
    self.rebuild_constraints_idx()
  
  def compress(self):
    return self.__getstate__()
  
  @classmethod
  def uncompress(cls,state):
    obj = cls.__new__(cls)
    obj.__setstate__(state)
    return obj
  
  def copy(self):
    copy_of_self = type(self).__new__(type(self))
    copy_of_self._nested_set_reprs = [r for r in self._nested_set_reprs]
    copy_of_self._distances_to_root = {k:v for k,v in self._distances_to_root.iteritems()}
    copy_of_self.free_leaves = {fl for fl in self.free_leaves}
    copy_of_self.score = self.score
    copy_of_self._pairs_accounted_for = {p for p in self._pairs_accounted_for}
    copy_of_self.built_clades = [bc for bc in self.built_clades]
    copy_of_self.constraints_idx = [c for c in self.constraints_idx]
    return copy_of_self
  
  def recompute(self,*args,**kwargs):
    if 'extension' in kwargs:
      extension = kwargs['extension']
      if type(extension).__name__ == 'LeafPairDistanceFrequency':
        for leaf in extension.leaves:
          self._distances_to_root[leaf] = 1
        self._pairs_accounted_for.add(extension.leaves)
        self._nested_set_reprs.append(frozenset({frozenset(extension.leaves),'r'}))
      else:
        self._distances_to_root = extension.distances_to_root
        self._pairs_accounted_for = extension.pairs_accounted_for
        self._nested_set_reprs = extension.nested_set_reprs
    else:
      if '_distances_to_root' in args:
        self._distances_to_root = {leaf:clade.trace_dist(leaf) for clade in self.built_clades
                                   for leaf in clade.leaf_names}
      if '_pairs_accounted_for' in args:
        self._pairs_accounted_for = {frozenset(pair) for clade in self.built_clades
                                     for pair in itertools.combinations(clade.leaf_names,2)}
      if '_nested_set_reprs' in args:
        self._nested_set_reprs = [frozenset({c.nested_set_repr(),'r'}) for c in self.built_clades]
  
  def _property_getter(self,property):
    try:
      return getattr(self,property)
    except AttributeError:
      self.recompute(property)
      return getattr(self,property)
  
  @property
  def current_clades_as_nested_sets(self):
    return self._property_getter('_nested_set_reprs')
  
  @property
  def distances_to_root(self):
    return self._property_getter('_distances_to_root')
  
  @property
  def pairs_accounted_for(self):
    return self._property_getter('_pairs_accounted_for')
  
  @property
  def complete(self):
    return len(self.built_clades) == 1 and not self.free_leaves
  
  def verify_remaining_proposed_pairs(self,extensions):
    for key,ext in extensions.items():
      for pair,dist in ext.unverified.items():
        # Looking up pair histogram separately in case lookup throws KeyError,
        # since we don't want to catch that one
        pair_histogram = self.pwdist_histograms_dict[pair]
        try:
          pair_freq = pair_histogram[dist]
        except KeyError:
          # If the resulting pairdist is not in the histogram, it means it wasn't observed
          # at all, so it's frequency is 0.
          pair_freq = 0.0
        if pair_freq < self.abs_cutoff:
          extensions.pop(key)
          break
        else:
          ext.score += math.log(pair_freq)
          ext.unverified.pop(pair)
          ext.verified.add(pair)
      else:
        # If all unverified pairs check out, make sure this extension is not passing
        # this filter entirely on the strength of pairs we just verified
        if not ext.consistent:
          extensions.pop(key)
    return extensions
  
  def as_nested_sets(self,extension):
    # Clades are represented by {clade.nested_set_repr(),'r'} (r for root) to indicate
    # their free-standing nature. This way {{clade1,'r'},{clade2,'r'}} represents two
    # free-standing clades, distinct from {clade1,clade2}, representing a single
    # free-standing clade (or tree) with two sub-clades.
    if hasattr(extension,'freq'):
      new_clade = frozenset({frozenset(extension.leaves),'r'})
      indeces_to_skip = []
    else:
      if hasattr(extension,'built_clade'):
        new_clade = frozenset({frozenset({extension.built_clade.clade.nested_set_repr(),
                                          extension.new_leaf}),'r'})
        indeces_to_skip = {extension.built_clade.index}
      else:
        new_clade = frozenset({frozenset(c.clade.nested_set_repr()
                                         for c in extension.clades),
                               'r'})
        indeces_to_skip = {c.index for c in extension.clades}
      extension.nested_set_reprs = [c for i,c in enumerate(self.current_clades_as_nested_sets)
                                    if i not in indeces_to_skip]
      extension.nested_set_reprs.append(new_clade)
    
    clades = [c for i,c in enumerate(self.current_clades_as_nested_sets)
                           if i not in indeces_to_skip]
    clades.append(new_clade)
    return clades
  
  @property
  def best_case(self):
    if not hasattr(self,'_best_case') or self._best_case is None:
      self._best_case = self.calculate_best_case()
    return self._best_case
  
  @property
  def nodes_left_to_build(self):
    if not hasattr(self,'_nodes_left_to_build')\
                                          or self._nodes_left_to_build is None:
      self._nodes_left_to_build = len(self.built_clades) +\
                                                      len(self.free_leaves) - 1
    return self._nodes_left_to_build
  
  @property
  def built_nodes_count(self):
    if not hasattr(self,'_built_nodes_count')\
                                            or self._built_nodes_count is None:
      self._built_nodes_count = self.total_nodes_to_build\
                                                     - self.nodes_left_to_build
    return self._built_nodes_count
  
  def reset(self):
    self._best_case = None
    self._nodes_left_to_build = None
    self._built_nodes_count = None
  
  @property
  def sort_key(self):
    return self.best_possible + self.score/len(self.pairs_accounted_for) if\
      float(self.built_nodes_count)/self.total_nodes_to_build < 0.4 else\
      self.best_case/self.built_nodes_count
  
  def calculate_best_case(self,pairs_accounted_for=None,distances_to_root=None,
                     score=None):
    pairs_accounted_for = pairs_accounted_for or self.pairs_accounted_for
    distances_to_root = distances_to_root or self.distances_to_root
    best_possible_final_score = score or self.score
    for pair,histogram in self.pwdist_histograms_dict.items():
      if pair not in pairs_accounted_for:
        min_dist = sum(distances_to_root[leaf] if leaf in distances_to_root else
                       0 for leaf in pair)+1
        acceptable_dist_freqs = [freq for dist,freq in histogram.items()
                                 if dist >= min_dist]
        if acceptable_dist_freqs:
          best_possible_final_score += math.log(max(acceptable_dist_freqs))
        else:
          return None
    return best_possible_final_score
  
  def best_case_with_extension(self,extension):
    try:
      updated_distances_to_root = {leaf:(self.distances_to_root[leaf]+1 if leaf in
                                         self.distances_to_root else 1)
                               for pair in itertools.chain(extension.consistent,
                                                           extension.verified)
                               for leaf in pair}
      for leaf in self.distances_to_root:
        if leaf not in updated_distances_to_root:
          updated_distances_to_root[leaf] = self.distances_to_root[leaf]
      
      updated_pairs_accounted_for = {pair for pair in
                                     itertools.chain(self.pairs_accounted_for,
                                                     extension.verified)}
      extension.distances_to_root = updated_distances_to_root
      extension.pairs_accounted_for = updated_pairs_accounted_for
    except AttributeError:
      updated_distances_to_root = dict(self.distances_to_root.iteritems())
      for leaf in extension.leaves:
        updated_distances_to_root[leaf] = 1
      updated_pairs_accounted_for = {pair for pair in self.pairs_accounted_for}
      updated_pairs_accounted_for.add(extension.leaves)
    
    try:
      best_possible_final_score = self.score + extension.score
    except AttributeError:
      best_possible_final_score = self.score + math.log(extension.freq)
    return self.calculate_best_case(updated_pairs_accounted_for,
                                    updated_distances_to_root,
                                    best_possible_final_score)
  
  def filter_proposed_extensions(self,new_pairs,joins,attachments,
                                 encountered,min_score=None):
    joins = self.verify_remaining_proposed_pairs(joins)
    attachments = self.verify_remaining_proposed_pairs(attachments)
    
    for extension_set in (new_pairs,joins,attachments):
      for key,extension in extension_set.items():
        nested_repr = self.as_nested_sets(extension)
        # First filter: has extension been encountered before?
        if encountered.already_encountered(nested_repr):
          extension_set.pop(key)
          continue
        # Second filter: is score with extension already worse than min_score?
        if min_score is not None:
          try: # Will fail if item is a new pair - forgiveness faster than permission
            if extension.score + self.score < min_score:
              extension_set.pop(key)
              continue
          except AttributeError: # Must be a new_pair:
            if math.log(extension.freq) + self.score < min_score:
              extension_set.pop(key)
              continue
        # Third filter: is there a way to extend the extension all the way to a full assembly?
        # Is the upper limit on best score for that assembly already worse than min_score?
        best_case = self.best_case_with_extension(extension)
        if best_case is None or best_case < min_score:
          extension_set.pop(key)
          continue
        encountered.remember(nested_repr)
    
    return new_pairs,joins,attachments
  
  def find_extensions(self,encountered,min_score=None):
    new_pairs = {}
    joins = self.KeyPassingDefaultDict(lambda key: ProposedExtension(*key))
    attachments = self.KeyPassingDefaultDict(lambda key: ProposedExtension(*key))
    already_connected = {frozenset(clade.leaf_names):i for i,clade in
                                                       enumerate(self.built_clades)}
    # All pairwise intersections should be empty
    assert not any(frozenset.intersection(*leafsetpair)
                   for leafsetpair in itertools.combinations(already_connected.iterkeys(),2))
    already_connected_splat = frozenset({leaf for leafset in already_connected
                                              for leaf in leafset})
    ac_leafdict = {leaf:self.IndexedClade(i,self.built_clades[i])
                   for leafset,i in already_connected.iteritems() for leaf in leafset}
    for i in self.constraints_idx:
      pair = self.constraints_master[i]
      if pair.dist == 1:
        # Pairs with distance 1 are added w/o questions. If continue with this path,
        # later we will make sure to remove from consideration all pairs that conflict this.
        new_pairs[i] = pair
      elif not pair.leaves & already_connected_splat:
        # If a pair has distance > 1 and neither leaf in pair has already been added to a
        # clade, then we can't do anything with it, so we silently skip it
        continue
      else:
        if pair.leaves <= already_connected_splat: # If both leaves are already in built clades
          if any(pair.leaves <= built_leafset for built_leafset in already_connected.iterkeys()):
            # If one built clade already contains both leaves, there is nothing to do either
            continue
          else:
            # If separate built clades contain one leaf each, then the pair goes into
            # the corresponding join's ProposedExtension object
            joins[frozenset(ac_leafdict[leaf] for leaf in pair.leaves)].check_pair(pair,i)
        else:
          # (pair.leaves & already_connected_splat) and (not pair.leaves <= already_connected_splat)
          # implies one leaf in pair is already contained in a built clade, the other isn't
          for leaf in pair.leaves: # Identify the new leaf and the built clade
            try:
              clade_of_attached_leaf = ac_leafdict[leaf]
              attached_leaf = leaf
            except KeyError as e:
              new_leaf = leaf
          # And put the pair into the corresponding attachment's ProposedExtension object
          attachments[frozenset({clade_of_attached_leaf,new_leaf})].check_pair(pair,i)
    return self.filter_proposed_extensions(new_pairs,joins,attachments,encountered,min_score)
    
  def build_extensions(self,new_pairs,joins,attachments):
    # Will need key of pair in new pairs, but not of keys in joins or attachments
    all_ext_to_build = new_pairs.items()+joins.values()+attachments.values()
    extended_assemblies = []
    while all_ext_to_build:
      extension = all_ext_to_build.pop()
      try:
        if not all_ext_to_build:
          extended_assemblies.append(extension.build_extension(self,in_place=True))
        else:
          extended_assemblies.append(extension.build_extension(self))
      except AttributeError:
        if not all_ext_to_build:
          build_in = self
        else:
          build_in = self.copy()
        idx_of_pair,pair = extension # Now we can get the key (index of pair in constraints_idx)
        
        # Select for dropping all pairs with distance 1 and one member of pair - they can't
        # have distance 1 with anyone except each other
        drop_these = [i for i in self.constraints_idx if self.constraints_master[i].dist == 1 and
                                              self.constraints_master[i].leaves & pair.leaves and
                                            not self.constraints_master[i].leaves == pair.leaves]
        # Select for dropping all pairs of these two leaves with distance > 1
        drop_these.extend(i for i in self.constraints_idx if self.constraints_master[i].dist > 1
                                            and self.constraints_master[i].leaves == pair.leaves)
        drop_these.append(idx_of_pair) # Finally, select for dropping this pair
        # Drop selected pairs from constraints_idx
        build_in.constraints_idx = filter(lambda x: x not in drop_these,build_in.constraints_idx)
        
        # Remove leaves in this pair from free_leaves
        for leaf in pair.leaves:
          build_in.free_leaves.remove(leaf) # Pop each leaf in pair from free_leaves
        
        # Build new clade and update the score
        build_in.built_clades.append(T.requisition(*tuple(pair.leaves)))
        build_in.recompute(extension=pair)
        build_in.score += math.log(pair.freq)
        extended_assemblies.append(build_in)
    for a in extended_assemblies:
      a.reset()
    return extended_assemblies
  
  def generate_extensions(self,encountered_assemblies,min_score=None):
    new_pairs,joins,attachments = self.find_extensions(encountered_assemblies,min_score)
    if any((new_pairs,joins,attachments)):
      return self.build_extensions(new_pairs, joins, attachments)
    else:
      return None


#===============================================================================
# Topology enumeration code
#===============================================================================

#------------------------------------------------------------------------------ 
# Global tracking of encoutered topology assembly states


class CladeReprTracker(object):
  def __init__(self,leaves):
    self.encountered = set()
    self.leaves = {leaf:i+1 for i,leaf in enumerate(sorted(leaves))}

  def __len__(self):
    return len(self.encountered)
  
  def _recursively_build_repr(self,c):
    nonleaflist = sorted((self._recursively_build_repr(m) for m in c if not
                          isinstance(m,str)),key=lambda x: x[1])
    nonleaves = ','.join(v[0] for v in nonleaflist)
    leafnumbers = sorted(self.leaves[m] for m in c if isinstance(m,str))
    leaves = ','.join(str(leaf) for leaf in leafnumbers)
    returnstr = '('
    take_my_min = []
    if nonleaves:
        returnstr += nonleaves
        take_my_min.append(nonleaflist[0][1])
        if leaves:
            returnstr += ','
    if leaves:
        returnstr += leaves
        take_my_min.append(leafnumbers[0])
    returnstr += ')'
    return returnstr,min(take_my_min)

  def make_str_repr(self,cladeset):
    cladestrlist = sorted((self._recursively_build_repr(m) for c in cladeset for m in c
                           if m !='r'),key=lambda x: x[1])
    return '['+','.join(v[0] for v in cladestrlist)+']'
  
  def already_encountered(self,cladeset):
    csrepr = self.make_str_repr(cladeset)
    if csrepr in self.encountered:
      return True
    else:
      self.encountered.add(csrepr)
      return False


class SharedCladeReprTracker(CladeReprTracker):
  def __init__(self,leaves,shared_dict):
    self.encountered = shared_dict
    self.leaves = {leaf:i+1 for i,leaf in enumerate(sorted(leaves))}
  
  def already_encountered(self, cladeset):
    return self.make_str_repr(cladeset) in self.encountered
  
  def remember(self,cladeset):
    self.encountered[self.make_str_repr(cladeset)] = None
  
  def forget(self,cladeset):
    try:
      self.encountered.pop(self.make_str_repr(cladeset))
    except KeyError:
      pass


#------------------------------------------------------------------------------ 
# Classes representing workspace in which topology construction takes place


class AssemblyWorkspace(object):
  def __init__(self,seed_assembly,num_requested_trees,max_workspace_size,
                    encountered_assemblies_storage,fifo=None,
                    track_min_score=True,acceptance_ratio_param=2.0,
                    acceptance_stiffness_param=1.0):
    if isinstance(seed_assembly,list):
      seed_assembly = TreeAssembly(*seed_assembly)
    self.workspace = [seed_assembly]
    self.accepted_assemblies = []
    self.rejected_assemblies = []
    self.encountered_assemblies = encountered_assemblies_storage
    
    self.num_requested_trees = num_requested_trees
    self._reached_num_requested_trees = False
    if track_min_score:
      self.curr_min_score = None
    self.max_workspace_size = max_workspace_size
    self.current_max = 10
    self.new_assembly_cache = []
    self.fifo = fifo
    
    self.push_cache = []
    self.iternum = 1
    self.total_nodes_to_build = seed_assembly.total_nodes_to_build
    self.push_count = 1 # "Pseudocount" to avoid division by 0
    self.topoff_count = 0
    self.topoff_param1 = 1
    self.topoff_param2 = 1
    self.accrp = acceptance_ratio_param
    self.accsp = acceptance_stiffness_param
  
  def check_if_num_requested_trees_reached(self):
    return len(self.accepted_assemblies) >= self.num_requested_trees
  
  @property
  def reached_num_requested_trees(self):
    if self._reached_num_requested_trees:
      return True
    elif self.check_if_num_requested_trees_reached():
      self._reached_num_requested_trees = True
      return True
    else:
      return False
  
  @property
  def acceptance_criterion(self):
    ratio = float(self.topoff_count)/self.push_count
    if ratio > self.accrp:
      return self.total_nodes_to_build
    elif ratio < 0.1:
      return 3
    else:
      return self.total_nodes_to_build - (self.total_nodes_to_build-3)*\
                          ((self.accrp-ratio)/(self.accrp-0.1))**self.accsp
  
  def apply_acceptance_logic_to_popped_assembly(self,popped,
                                                     rejected_assemblies,
                                                     counter,no_reject=False):
    if popped[2] > self.curr_min_score:
      self.topoff_count += 1
      if popped[3] <= self.acceptance_criterion or no_reject == True:
        uncompressed_assembly = TreeAssembly.uncompress(popped)
        self.log("TopoffAccepted",uncompressed_assembly,best_case=popped[2])
        self.workspace.append(uncompressed_assembly)
      else:
        self.log("TopoffPostponed",popped,compressed=True)
        counter[0] += 1
        rejected_assemblies.append(popped)
    else:
      uncompressed_assembly = TreeAssembly.uncompress(popped)
      self.encountered_assemblies.forget(
                           uncompressed_assembly.current_clades_as_nested_sets)
      self.log("TopoffRejected",uncompressed_assembly)

  
  def fill_workspace_from_fifo(self,max_size,rejected_assemblies,counter):
    while len(self.workspace) < max_size and counter[0] < 100:
      popped = self.fifo.pop()
      if popped is None:
        break
      else:
        self.apply_acceptance_logic_to_popped_assembly(popped,
                                                       rejected_assemblies,
                                                       counter)
  
  def top_off_workspace(self,max_size=None):
    if self.fifo is None or isinstance(self.fifo,str):
      return
    max_size = max_size or self.max_workspace_size
    while len(self.workspace) < max_size:
      counter = [0]
      rejected = []
      self.fill_workspace_from_fifo(max_size,rejected,counter)
      self.push_to_fifo(rejected)
      if counter[0] < 100:
        break
  
  def push_to_fifo(self,push_these):
    if self.fifo is None:
      self.fifo = FIFOfile()
    elif isinstance(self.fifo,str):
      self.fifo = FIFOfile(name=self.fifo)
    for item in push_these:
      self.fifo.push(item)
  
  def purge_push_cache(self):
    self.push_to_fifo(self.push_cache)
    self.push_cache = []
  
  def push(self,*args):
    self.push_cache.extend(item.compress() for item in args)
    for assembly in args:
      self.log("CachedToPush",assembly)
    if len(self.push_cache) > 100:
      self.purge_push_cache()
    self.push_count += len(args)
  
  def sock_away_extras(self,too_many_here,max_size=None):
    max_size = max_size or self.max_workspace_size
    while len(too_many_here) > max_size:
      self.push(too_many_here.pop())
  
  def update_workspace(self,new_assemblies):
#     for assembly in new_assemblies:
#       self.log("CachedNewAssembly",assembly)
    self.new_assembly_cache.extend(new_assemblies)
    self.new_assembly_cache.sort(key=lambda a:a.sort_key,reverse=True)
    max_size = self.max_workspace_size if self.reached_num_requested_trees\
                                                          else self.current_max
    if len(self.new_assembly_cache) > max_size:
      self.sock_away_extras(self.new_assembly_cache,max_size)
  
  def finalize_workspace(self):
    self.workspace.extend(self.new_assembly_cache)
    self.workspace.sort(key=lambda a: a.sort_key,reverse=True)
    self.new_assembly_cache = []
    if not self.reached_num_requested_trees:
      if self.workspace:
        small_wksp_criterion = max(a.nodes_left_to_build
                                   for a in self.workspace) > 3
      else:
        small_wksp_criterion = True
      self.current_max = 10 if small_wksp_criterion\
                                            else min(self.max_workspace_size,
                                                     100)
      if len(self.workspace) > self.current_max:
        self.sock_away_extras(self.workspace,self.current_max)
      if len(self.workspace) < self.max_workspace_size:
        self.top_off_workspace(self.current_max)
    else:
      self.current_max = min(self.max_workspace_size,max(10,
                             self.max_workspace_size/(50*max(0.02,
                                                           (self.topoff_param1/
                                                            self.topoff_param2)
                                                                           ))))
      if self.workspace:
        if len(self.workspace) > self.max_workspace_size:
          self.sock_away_extras(self.workspace)
        if len(self.workspace) < self.max_workspace_size:
          self.top_off_workspace()
        self.topoff_param1 += 1.0
        self.topoff_param2 = max(1.0,self.topoff_param2 - 1.0)
      else:
        self.top_off_workspace(self.current_max)
        self.topoff_param1 = 1.0
        self.topoff_param2 += 5.0
    self.purge_push_cache()
  
  def check_completion_status(self,assembly):
    if assembly.complete:
      if not self.reached_num_requested_trees\
                                       or assembly.score > self.curr_min_score:
        self.accepted_assemblies.append(assembly)
        self.accepted_assemblies.sort(key=lambda x: x.score,reverse=True)
        while len(self.accepted_assemblies) > self.num_requested_trees:
          self.rejected_assemblies.append(self.accepted_assemblies.pop())
        self.curr_min_score = self.accepted_assemblies[-1].score
      else:
        self.log("CompleteRejected",assembly)
        self.rejected_assemblies.append(assembly)
      return
    else:
      return assembly
  
  def prepare_to_terminate(self):
    self.workspace.extend(self.new_assembly_cache)
    while self.workspace:
      self.push(self.workspace.pop())
      self.purge_push_cache()
    self.ready_to_terminate = True
  
  def iterate(self,interrupt_callable=lambda: False):
    drop_from_workspace_idx = []
    workspace_this_iter = [assembly for assembly in self.workspace]
    for i,assembly in enumerate(workspace_this_iter):
      if interrupt_callable():
        continue
      self.log("WorkingOn",assembly)
      if assembly.best_case < self.curr_min_score:
        self.log("AbandonedBestCase",assembly)
        drop_from_workspace_idx.append(i)
        continue
      else:
        extended_assemblies = assembly.generate_extensions(
                                                   self.encountered_assemblies,
                                                           self.curr_min_score)
        if extended_assemblies is None:
          self.log("AbandonedNoExtensions",assembly)
          drop_from_workspace_idx.append(i)
          continue
        else:
          assert assembly is extended_assemblies[-1]
          if self.check_completion_status(extended_assemblies.pop()) is None:
            drop_from_workspace_idx.append(i)
        self.update_workspace([asbly for asbly in extended_assemblies
                               if self.check_completion_status(asbly)
                                                                  is not None])
    for i in drop_from_workspace_idx[::-1]:
      if self.workspace[i] in self.accepted_assemblies:
        # Don't want to forget an accepted assembly because that can result in
        # accepting the same assembly twice!
        self.workspace.pop(i)
      else:
        self.encountered_assemblies.forget(
                           self.workspace.pop(i).current_clades_as_nested_sets)
    if interrupt_callable():
      self.prepare_to_terminate()
    else:
      self.finalize_workspace()
    self.iternum += 1


class WorkerProcAssemblyWorkspace(AssemblyWorkspace):
  class AssemblyWorkFinished(Exception):
    pass
  
  def __init__(self,fifo,queue,min_score,shared_encountered_assemblies_dict,
               score_submission_queue,start_time_val,leaves_to_assemble,seed_assembly,
               num_requested_trees,max_workspace_size,
               max_monitor_file_size=100*1024**2,**kwargs):
    encountered_assemblies = SharedCladeReprTracker(leaves_to_assemble,
                                            shared_encountered_assemblies_dict)
    AssemblyWorkspace.__init__(self,seed_assembly,num_requested_trees,
                                    max_workspace_size,encountered_assemblies,
                                    fifo,track_min_score=False,**kwargs)
    
    self._curr_min_score = min_score
    self.queue = queue
    self.score_submission_queue = score_submission_queue
    
    self.start_time = start_time_val
    self._max_monitor_file_size = max_monitor_file_size
    self._monitor_file_count = 1
    self.proc_name = multiprocessing.current_process().name
  
  def check_if_num_requested_trees_reached(self):
    # Multiply initial value by 0.9, because who knows if the comparison
    # -sys.float_info.max > -sys.float_info.max may sometimes return true?
    return self._curr_min_score.value > -sys.float_info.max*0.9

  @property
  def monitor(self):
    if not hasattr(self,'_monitor'):
      self._monitor = open(self.proc_name+'_activity_dump001','w',0)
    elif os.path.getsize(self._monitor.name) > self._max_monitor_file_size:
      self._monitor.close()
      self._monitor_file_count += 1
      self._monitor = open(self.proc_name+'_activity_dump'+\
                           str(self._monitor_file_count).zfill(3),
                           'w',0)
    return self._monitor
  
  @property
  def complete_trees_fh(self):
    if not hasattr(self,'_complete_trees_fh'):
      self._complete_trees_fh = open(multiprocessing.current_process().name+\
                                     '_complete_trees','w',0)
    return self._complete_trees_fh

  @property
  def curr_min_score(self):
    return self._curr_min_score.value
  
  @property
  def time_stamp(self):
    return "%0.5f" % (time.time()-self.start_time.value)
  
  def log(self,message,assembly=None,proc_stamp=True,compressed=False,best_case=None):
    currscore = self.curr_min_score
    currscore = 'min_score='+repr(None) if currscore == -sys.float_info.max\
                                           else "min_score=%0.5f" % currscore
    iternum = "iter="+str(self.iternum)
    pushcount = "push_count="+str(self.push_count)
    topoffcount = "topoff_count="+str(self.topoff_count)
    criterion = "criterion=%0.2f" % self.acceptance_criterion
    if proc_stamp:
      stamp = 'STAMP: '+'  '.join([multiprocessing.current_process().name,
                                  'time='+self.time_stamp,currscore,iternum,
                                  pushcount,topoffcount,criterion])+'   '
    else:
      stamp = 'STAMP: '+'  '.join(['time='+self.time_stamp,currscore,iternum,
                                  pushcount,topoffcount,criterion])+'   '
    if assembly is None:
      print >>self.monitor,stamp,message
    elif compressed:
      print >>self.monitor,stamp,"\t%0.5f\t" % assembly[2],assembly[3],\
                                 "\t%0.5f\t" % (assembly[2]/assembly[3]),\
                                 message
    else:
      best_case = best_case or assembly.best_case
      print >>self.monitor,stamp,"\t%0.5f\t" % best_case,\
                                 assembly.nodes_left_to_build,\
                                 "\t%0.5f\t" % assembly.sort_key,message,\
                                 self.encountered_assemblies.make_str_repr(
                                       assembly.current_clades_as_nested_sets)
  
  def fill_workspace_from_fifo(self,max_size,rejected_assemblies,counter):
    no_reject = False
    while len(self.workspace) < max_size and counter[0] < 100:
      if counter[0] == 99 and self.topoff_count == counter[0]:
        counter[0] = 0
        no_reject = True
      try:
        pickled_assembly = self.queue.get_nowait()
      except Queue.Empty:
        if self.workspace:
          break
        if rejected_assemblies:
          while rejected_assemblies:
            self.push_to_fifo((rejected_assemblies.pop(),))
          self.purge_push_cache()
          counter[0] = 0
          no_reject = True
          time.sleep(5)
        else:
          try:
            pickled_assembly = self.queue.get(timeout=5)
          except Queue.Empty:
            if not self.fifo.is_set():
              raise self.AssemblyWorkFinished
            else:
              continue
      self.apply_acceptance_logic_to_popped_assembly(pickled_assembly,
                                                     rejected_assemblies,
                                                     counter,no_reject)
  
  def push_to_fifo(self,push_these):
#     for item in push_these:
#       self.log("Pushing",item,compressed=True)
    self.fifo.push_all(item for item in push_these)
  
  def check_completion_status(self,assembly):
    if assembly.complete:
      if assembly.score > self.curr_min_score:
        self.log("CompleteAccepted",assembly)
        self.score_submission_queue.put(assembly.score)
        self.accepted_assemblies.append(assembly)
        self.complete_trees_fh.write(str(assembly.score)+'\t'+assembly.built_clades[0].write('as_string','newick',plain=True))
      else:
        self.log("CompleteRejected",assembly)
        self.rejected_assemblies.append(assembly)
      for i in xrange(len(self.accepted_assemblies)-1,-1,-1):
        if self.accepted_assemblies[i].score < self.curr_min_score:
          self.rejected_assemblies.append(self.accepted_assemblies.pop(i))
      return
    else:
      self.log("Extended",assembly)
      return assembly
  
  def iterate(self,*args,**kwargs):
    try:
      print >>self.monitor,'-'*80
      print >>self.monitor,"START OF ITERATION",self.iternum,
      print >>self.monitor,"\tworkspace size:",len(self.workspace),
      print >>self.monitor,"\ttime:",self.time_stamp
      AssemblyWorkspace.iterate(self,*args,**kwargs)
    except self.AssemblyWorkFinished:
      print >>self.monitor,"Caught 'FINISHED' signal"
      self.iternum += 1
      return 'FINISHED'
    finally:
      print >>self.monitor,"END OF ITERATION",self.iternum-1,
      print >>self.monitor,"\tworkspace size:",len(self.workspace),
      print >>self.monitor,"\ttime:",self.time_stamp
      print >>self.monitor,'-'*80
      gc.collect()


#------------------------------------------------------------------------------
# Process classes which load incomplete assemblies into queue


class QueueLoader(multiprocessing.Process):
  def __init__(self,fifo,close_fifo_EV,queue,interrupt,start_loading):
    multiprocessing.Process.__init__(self,name=multiprocessing.current_process().name+'--QueueLoader')
    self.fifo = fifo
    self.close_fifo = close_fifo_EV
    self.queue = queue
    self.interrupt = interrupt
    self.start_loading = start_loading
    self.daemon = True
  
  def put_in_queue(self,item):
    while True:
      try:
        self.queue.put(item,timeout=5)
        break
      except Queue.Full:
        continue
  
  def run(self):
    self.fifo.start_OUT_end()
    self.start_loading.wait()
    while not self.close_fifo.is_set():
      popped = self.fifo.pop()
      if popped is None:
        if self.interrupt.is_set():
          self.put_in_queue('FIFO_EMPTY')
          break
        else:
          continue
      else:
        self.put_in_queue(popped)
    self.fifo.close()
    return


def RestartQueueReloader_worker_task(stored_string):
  return RestartQueueReloader.prepare_assembly_state(stored_string,
                                                   globals()['leaf_name_map'],
                                                   globals()['dummy_assembly'])
  
class RestartQueueReloader(multiprocessing.Process):
  def __init__(self,restart_from,queue,release_loaders,dummy_assembly,
                    initial_batch_size=1000,numproc=None):
    multiprocessing.Process.__init__(self,name='RestartQueueLoader')
    self.restart_from = restart_from
    self.queue = queue
    self.start_workers = multiprocessing.Event()
    self.release_loaders = release_loaders
    self.initial_batch_size = initial_batch_size
    self.dummy_assembly = dummy_assembly
    self.leaf_name_map = {i+1:leaf for i,leaf in enumerate(sorted(
                                                dummy_assembly.leaves_master))}
    self.numproc = numproc
  
  @staticmethod
  def pool_worker_init(leaf_name_map,dummy_assembly):
    globals()['leaf_name_map'] = leaf_name_map
    globals()['dummy_assembly'] = dummy_assembly
  
  @staticmethod
  def prepare_assembly_state(stored_string,leaf_name_map,dummy_assembly):
    built_clades,score,best_case,left_to_build = eval(stored_string)
    built_clades = [T_BASE(StringIO(c)) for c in built_clades]
    for c in built_clades:
      for l in c.get_terminals():
        l.name = leaf_name_map[eval(l.name)]
    built_repr = ';'.join(''.join(repr(dummy_assembly.convert_containers(
                                                 c.nested_set_repr())).split())
             for c in built_clades)
    for k,v in dummy_assembly.pickle_encoding.items():
      built_repr = built_repr.replace(v,k)
    return built_repr,score,best_case,left_to_build
  
  def run(self):
    with tarfile.open(self.restart_from) as tf:
      fh = tf.extractfile('./unfinished_assemblies')
      worker_pool = multiprocessing.Pool(self.numproc,
                                         initializer=self.pool_worker_init,
                                         initargs=(self.leaf_name_map,
                                                   self.dummy_assembly))
      decoded_assemblies = worker_pool.imap(RestartQueueReloader_worker_task,
                                            fh)
      worker_pool.close()
      for i,prepared_state in enumerate(decoded_assemblies):
        try:
          self.queue.put(prepared_state,timeout=5)
        except Queue.Full:
          break
        if i+1 == self.initial_batch_size:
          break
      else:
        self.start_workers.set()
        self.release_loaders.set()
        fh.close()
        return
      self.start_workers.set()
      for prepared_state in decoded_assemblies:
        while True:
          try:
            self.queue.put(prepared_state,timeout=5)
            break
          except Queue.Full:
            continue
      self.release_loaders.set()
      fh.close()
      worker_pool.join()
      return


#------------------------------------------------------------------------------ 
# Process class in which parallelized work takes place


class AcceptedAssembly(object):
  __slots__ = ('score', 'topology')
  
  def __init__(self,score=None,topology=None,string_repr=None):
    if string_repr is not None:
      score_str,topology = string_repr.split()
      self.score = float(score_str)
      self.topology = topology
    else:
      self.score = score
      self.topology = topology.strip()
  
  def __repr__(self):
    return str(self.score)+'\t'+self.topology+'\n'

class AssemblerProcess(multiprocessing.Process):
  def __new__(cls,*args,**kwargs):
    if not hasattr(cls,'instcount'):
      cls.instcount = 0
    cls.instcount += 1
    return multiprocessing.Process.__new__(cls,*args,**kwargs)
  
  def __init__(self,queue,shared_encountered_assemblies_dict,shared_min_score,
                    score_submission_queue,seed_assembly,pass_to_workspace,
                    start_time_val,results_queue,release_queue_loader,
                    fifo_max_file_size=1.0):
    multiprocessing.Process.__init__(self,name='AssemblerProcess-'\
                                                 +str(self.instcount).zfill(3))
    
    self.queue = queue
    self.min_score = shared_min_score
    self.encountered_assemblies_dict = shared_encountered_assemblies_dict
    self.score_submission_queue = score_submission_queue
    self.pass_to_workspace = pass_to_workspace
    self.seed_assembly = seed_assembly
    self.fifo_max_file_size = fifo_max_file_size
    
    self.start_time = start_time_val
    
    self.close_fifo = multiprocessing.Event()
    
    self.interrupt = multiprocessing.Event()
    self.QueueLoader_interrupt = multiprocessing.Event()
    self.finished = multiprocessing.Event()
    self.shutdown = multiprocessing.Event()
    self.release_queue_loader = release_queue_loader
    self.results_queue = results_queue
  
  def enqueue_results(self):
    for assembly in self.assemblies.accepted_assemblies:
      accepted = AcceptedAssembly(assembly.score,
                                  assembly.built_clades[0].write('as_string',
                                                               format='newick',
                                                                 plain=True))
      self.results_queue.put(accepted)
    self.results_queue.put('FINISHED')
  
  def run(self):
    self.fifo = SharedFIFOfile(suffix='--'+self.name,
                               max_file_size_GB=self.fifo_max_file_size)
    self.pass_to_workspace.kwargs['seed_assembly'] = self.seed_assembly
    self.assemblies = WorkerProcAssemblyWorkspace(self.fifo,self.queue,self.min_score,
                                                  self.encountered_assemblies_dict,
                                                  self.score_submission_queue,
                                                  self.start_time,
                                                  *self.pass_to_workspace.args,
                                                  **self.pass_to_workspace.kwargs)
    self.queue_loader_p = QueueLoader(self.fifo,self.close_fifo,self.queue,
                                      self.QueueLoader_interrupt,
                                      self.release_queue_loader)
    self.queue_loader_p.start()
    self.fifo.start_IN_end()
    try:
      iter_result = None
      iter_counter = 0
      while not self.shutdown.is_set():
        iter_counter += 1
        iter_result = self.assemblies.iterate(self.interrupt.is_set)
        if self.interrupt.is_set():
          if not hasattr(self.assemblies,'ready_to_terminate'):
            self.assemblies.prepare_to_terminate()
          self.QueueLoader_interrupt.set()
          self.fifo.set()
          self.enqueue_results()
          self.shutdown.wait()
          break
        if iter_result == 'FINISHED':
          self.finished.set()
          if self.shutdown.wait(5):
            self.enqueue_results()
            break
          else:
            continue
      self.assemblies.monitor.close()
      self.assemblies.complete_trees_fh.close()
      self.close_fifo.set()
      self.fifo.close()
      self.queue_loader_p.join(timeout=15)
    except:
      self.fifo.current_writing_file.wh.close()
      self.fifo.tmpdir_obj.__exit__(None,None,None)
      self.queue_loader_p.terminate()
      raise
    return


#------------------------------------------------------------------------------ 
# Central enumeration process class


def unpickle_assembly_for_save(pickled_assembly_state):
  state = globals()['zeroth_assembly']._unpack_state(pickled_assembly_state)
  state['built_clades'] = [T.rebuild_on_unpickle(c).write('as_string',
                                                           format='newick',
                                                           plain=True)
                           for c in state['built_clades']]
  for i,c in enumerate(state['built_clades']):
    for k,v in globals()['leaf_name_encoding']:
      c = c.replace(k,str(v))
    state['built_clades'][i] = c
  return state


class MainTopologyEnumerationProcess(multiprocessing.Process):
  def __init__(self,leafdist_histograms,constraint_freq_cutoff=0.9,
                    absolute_freq_cutoff=0.01,max_workspace_size=10000,
                    max_queue_size=10000,fifo_max_file_size=1.0,
                    num_requested_topologies=1000,num_workers=None,
                    save_file_name='early_termination_save',
                    restart_from=None,**kwargs):
    multiprocessing.Process.__init__(self)
    self.assembly_queue_manager = multiprocessing.Manager()
    self.assembly_queue = self.assembly_queue_manager.Queue(max_queue_size)
    self.encountered_assemblies_manager = multiprocessing.Manager()
    self.encountered_assemblies_dict = self.encountered_assemblies_manager.dict()
    self.results_queue = multiprocessing.Queue()
    self.scores_queue = multiprocessing.Queue()
    self.min_score = multiprocessing.Value('d',-sys.float_info.max)
    self.get_PIDs,self.send_PIDs = multiprocessing.Pipe(duplex=False)
    self.release_queue_loaders = multiprocessing.Event()
    self.stop = multiprocessing.Event()
    self.save_written = multiprocessing.Event()
    self.finished = multiprocessing.Event()
    self.shutdown = multiprocessing.Event()
    
    self.start_time = multiprocessing.Value('d',time.time())
    
    self.histograms = leafdist_histograms
    self.leaves = {l for pair in self.histograms for l in pair[0]}
    self.constraint_freq_cutoff = constraint_freq_cutoff
    self.absolute_freq_cutoff = absolute_freq_cutoff
    self.max_workspace_size = max_workspace_size
    self.fifo_max_file_size = fifo_max_file_size
    self.num_requested_topologies = num_requested_topologies
    self.num_workers = multiprocessing.cpu_count() if num_workers is None\
                                                               else num_workers
    self.expected_number_results_queue_sentinels = self.num_workers
    self.zeroth_assembly = TreeAssembly(self.histograms,
                                        self.constraint_freq_cutoff,
                                        self.leaves,self.absolute_freq_cutoff,
                                        keep_alive_when_pickling=False)
    self.save_file_name = save_file_name
    self.restart_from = restart_from
    if restart_from is not None:
      self.initial_batch_size = max(int(round(max_queue_size*0.1)),1000)
      self.expected_number_results_queue_sentinels += 1
    self.kwargs = kwargs
  
  def clean_up(self):
    self.assembly_queue_manager.shutdown()
    self.encountered_assemblies_manager.shutdown()
    self.results_queue.close()
    self.results_queue.join_thread()
    self.scores_queue.close()
    self.scores_queue.join_thread()
    self.send_PIDs.close()
    self.get_PIDs.close()
  
  def assemblies_from_queue_generator(self):
    emptied_FIFO_counter = 0
    while emptied_FIFO_counter < self.num_workers:
      try:
        pickled_assembly_state = self.assembly_queue.get(timeout=5)
        if pickled_assembly_state == 'FIFO_EMPTY':
          emptied_FIFO_counter += 1
        else:
          yield pickled_assembly_state
      except Queue.Empty:
        continue 
  
  @staticmethod
  def pool_assembly_unpickler_worker_init(zeroth_assembly,leaf_name_encoding):
    globals()['zeroth_assembly'] = zeroth_assembly
    globals()['leaf_name_encoding'] = leaf_name_encoding.items() 
  
  def write_save(self):
    leaf_name_encoding = CladeReprTracker(self.leaves).leaves
    reverse_encoding = {v:k for k,v in leaf_name_encoding.items()}
    os.mkdir('tmp_savedir')
    with open('tmp_savedir/leaf_name_encoding','w',0) as wh:
      wh.write(repr(reverse_encoding)+'\n')
    worker_pool = multiprocessing.Pool(self.num_workers,
                          initializer=self.pool_assembly_unpickler_worker_init,
                                       initargs=(self.zeroth_assembly,
                                                 leaf_name_encoding))
    assembly_states = worker_pool.imap(unpickle_assembly_for_save,
                                       self.assemblies_from_queue_generator())
    worker_pool.close()
    with open('tmp_savedir/unfinished_assemblies','w',0) as wh:
      for state in assembly_states:
        wh.write(repr((state['built_clades'],round(state['score'],5),
                       round(state['_best_case'],5),
                       state['_nodes_left_to_build'])).replace(' ','')+'\n')
    worker_pool.join()
    with open('tmp_savedir/encountered_assemblies','w',0) as wh:
      while True:
        try:
          assembly_repr,_ = self.encountered_assemblies_dict.popitem()
          wh.write(assembly_repr+'\n')
        except KeyError:
          break
    
    accepted = []
    finished_worker_counter = 0
    while finished_worker_counter < self.expected_number_results_queue_sentinels:
      received = self.results_queue.get()
      if received == 'FINISHED':
        finished_worker_counter += 1
      else:
        accepted.append(received)
    accepted.sort(key=lambda x: x.score,reverse=True)
    with open('tmp_savedir/accepted_complete_assemblies','w',0) as wh:
      for assembly in accepted:
        wh.write(repr(assembly))
    import shutil
    shutil.make_archive(self.save_file_name,'gztar','tmp_savedir')
    shutil.rmtree('tmp_savedir')
  
  def set_up_initial_run(self,workspace_args):
    self.release_queue_loaders.set()
    seed_assemblies = [self.zeroth_assembly]
    while len(seed_assemblies) < self.num_workers:
      seed_assemblies = [a for old_a in seed_assemblies
                         for a in old_a.generate_extensions(
                                  SharedCladeReprTracker(self.leaves,
                                            self.encountered_assemblies_dict))]
    seed_assemblies.sort(key=lambda a: a.sort_key)
    procs = [AssemblerProcess(self.assembly_queue,
                              self.encountered_assemblies_dict,
                              self.min_score,self.scores_queue,
                              seed_assemblies.pop(),workspace_args,
                              self.start_time,self.results_queue,
                              self.release_queue_loaders,
                              self.fifo_max_file_size)
             for i in xrange(self.num_workers)]
    while seed_assemblies:
      self.assembly_queue.put(seed_assemblies.pop().compress())
    self.accepted_scores = []
    return procs
  
  def set_up_restarted_run(self,workspace_args):
    self.restart_queue_loader = RestartQueueReloader(self.restart_from,
                                                     self.assembly_queue,
                                                     self.release_queue_loaders,
                                                     self.zeroth_assembly,
                                                     self.initial_batch_size,
                                                     self.num_workers)
    self.restart_queue_loader.start()
    with tarfile.open(self.restart_from) as tf:
      fh = tf.extractfile('./accepted_complete_assemblies')
      self.previously_accepted_assemblies = []
      for l in fh:
        self.previously_accepted_assemblies.append(
                                               AcceptedAssembly(string_repr=l))
      self.previously_accepted_assemblies.sort(key=lambda x: x.score,
                                               reverse=True)
      self.accepted_scores = [a.score
                                  for a in self.previously_accepted_assemblies]
      while len(self.accepted_scores) > self.num_requested_topologies:
        self.accepted_scores.pop()
        self.previously_accepted_assemblies.pop()
      if len(self.accepted_scores) == self.num_requested_topologies:
        self.min_score.value = self.accepted_scores[-1]
      fh.close()
      for assembly in self.previously_accepted_assemblies:
        self.results_queue.put(assembly)
      self.results_queue.put('FINISHED')
      fh = tf.extractfile('./leaf_name_encoding')
      reverse_leaf_map = {v:k for k,v in CladeReprTracker(self.leaves).leaves.items()}
      assert eval(fh.read()) == reverse_leaf_map
      fh.close()
      fh = tf.extractfile('./encountered_assemblies')
      for l in fh:
        self.encountered_assemblies_dict[l.strip()] = None
      fh.close()
    self.restart_queue_loader.start_workers.wait()
    procs = [AssemblerProcess(self.assembly_queue,
                              self.encountered_assemblies_dict,
                              self.min_score,self.scores_queue,
                              TreeAssembly.uncompress(self.assembly_queue.get()),
                              workspace_args,
                              self.start_time,self.results_queue,
                              self.release_queue_loaders,
                              self.fifo_max_file_size)
             for i in xrange(self.num_workers)]
    return procs
  
  def run(self):
    self.kwargs.update({'num_requested_trees':self.num_requested_topologies,
                        'max_workspace_size':self.max_workspace_size})
    workspace_args = namedtuple('ArgsKwargs',
                                ['args','kwargs'])((self.leaves,),self.kwargs)
    try:
      if self.restart_from is not None:
        self.procs = self.set_up_restarted_run(workspace_args)
      else:
        self.procs = self.set_up_initial_run(workspace_args)
      for p in self.procs:
        p.start()
        self.send_PIDs.send((p.pid,p.name))
      
      while any(not p.finished.is_set() for p in self.procs):
        if self.stop.is_set():
          for p in self.procs:
            p.interrupt.set()
          self.write_save()
          self.save_written.set()
          break
        try:
          proposed_score = self.scores_queue.get(timeout=0.05)
          if len(self.accepted_scores) < self.num_requested_topologies\
                                  or proposed_score > self.accepted_scores[-1]:
            self.accepted_scores.append(proposed_score)
            self.accepted_scores.sort(reverse=True)
            while len(self.accepted_scores) > self.num_requested_topologies:
              self.accepted_scores.pop()
            if len(self.accepted_scores) == self.num_requested_topologies:
              self.min_score.value = self.accepted_scores[-1]
        except Queue.Empty:
          continue
      
      for p in self.procs:
        p.shutdown.set()
      
      self.finished.set()
      while not self.shutdown.wait(1):
        continue
      for p in self.procs:
        p.join(timeout=10)
        if p.is_alive():
          p.terminate()
          p.join()
      if hasattr(self,'restart_queue_loader'):
        self.restart_queue_loader.join(timeout=10)
        if self.restart_queue_loader.is_alive():
          self.restart_queue_loader.terminate()
          self.restart_queue_loader.join()
    
    except Exception:
      if hasattr(self,'procs'):
        for p in self.procs:
          p.shutdown.set()
          p.join(timeout=15)
          if p.is_alive():
            p.terminate()
            p.join()
      raise

