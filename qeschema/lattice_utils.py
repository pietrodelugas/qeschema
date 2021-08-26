#
# Copyright (c), 2015-2021, Quantum Espresso Foundation and SISSA (Scuola
# Internazionale Superiore di Studi Avanzati). All rights reserved.
# This file is distributed under the terms of the MIT License. See the
# file 'LICENSE' in the root directory of the present distribution, or
# http://opensource.org/licenses/MIT.
#
# Authors: Davide Brunato, Pietro Delugas
#

import logging 
import numpy as np
from math import isclose

logger = logging.getLogger("qeschema")

class lattice:
  def __init__(self, vectors = None, abc = None, abc_is_primitive = False, bravais_index = None, alt_axes = None, celldm = None):
    self.cell  = vectors
    if abc_is_primitive:
      self.primitive_abc = abc
      self.conventional_abc = None
    else:
      self.conventional_abc = abc
      self.primitive_abc = None
    self.bravais_index = bravais_index
    self.alt_axes = alt_axes
    self.celldm = celldm 
    if vectors:
      self.init_from_vectors()
    elif abc:
      self.init_from_abc() 
    elif celldm:
      self.init_from_celldm() 
    
  def init_from_vectors(self):





def abc_from_cell(cell, bravais_index, alternative_axes = None, dgts=5):
  """
  cell: dictionary with 'a1','a2','a3' keys with values corresponding to the 3 lattice vectors
  returns: 6ple with values of A,B,C,cos(AB), cos(BC), cos(AC) lenghts in Angstrom and cosines of the 
           lattice vectors
  """
  at = np.array((cell['a1'],cell['a2'], cell['a3']))
  celldm, signed_ibrav = at2celldm(at, ibrav)
  at2   = latgen(signed_ibrav, celldm)
  check = compare_at(at, at2, dgts = dgts) 
  if check: 
    a,b,c,cosab,cosbc,cosca = celldm2abc(celldm, ibrav,dgts=dgts) 
  else:
    logger.error(f"Cell parameters in atomic_structure/cell are incompatible with bravais_index={bravais_index}")
  return ibrav, (a, b, c, cosab, cosbc, cosca) 

def at2celldm(a, ibrav):
  """
  a: numpy array with the cell parameters
  ibrav: bravais index 
  """ 
  c1 = 0.0
  c2 = 0.0 
  c3 = 0.0
  c4 = 0.0
  c5 = 0.0
  c6 = 0.0  
  if ibrav == 0:
    c1 = np.sqrt(a[0].dot(a[0]))
  elif ibrav == 1:
