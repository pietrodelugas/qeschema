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

def abc_from_cell(cell, bravais_index, dgts=5):
  """
  cell: dictionary with 'a1','a2','a3' keys with values corresponding to the 3 lattice vectors
  returns: 6ple with values of A,B,C,cos(AB), cos(BC), cos(AC) lenghts in Angstrom and cosines of the 
           lattice vectors
  """
  at = np.array((cell['a1'],cell['a2'], cell['a3']))
  celldm, ibrav = at2celldm(at,ibrav)
  at2   = latgen(ibrav, celldm)
  check = compare_at(at,at2,dgts=dgts) 
  if check: 
    a,b,c,cosab,cosbc,cosca = celldm2abc(celldm, ibrav,dgts=dgts) 
  else:
    logger.error(f"Cell parameters in atomic_structure/cell are incompatible with bravais_index={bravais_index}")
  return ibrav, (a, b, c, cosab, cosbc, cosca) 

def at2cell(a, ibrav):
  """
  a: numpy array with the cell parameters
  ibrav: bravais index 
  """ 
  c1 = np.sqrt(a[0].dot(a[0]))
  c2 = np.sqrt(a[1].dot(a[1]))/c1
  c3 = np.sqrt(a[2].dot(a[2]))/c1
  c4 = np.
  if ibrav == 0:
    
