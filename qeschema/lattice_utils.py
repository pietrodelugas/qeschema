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

from xmlschema.validators.helpers import base64_binary_validator

logger = logging.getLogger("qeschema")


def qe_ibrav(bravais_index, alt_axes = None):
  if bravais_index < 0 or bravais_index > 14:
    logger.error("wrong value for the bravais_index") 
  if bravais_index == 3:
    if alt_axes == "b:a-b+c:-c": 
      ibrav_ =  -3 
  elif bravais_index == 5:
    if alt_axes == "3fold-111":
      ibrav_ =  -5 
  elif bravais_index == 9:
    if alt_axes == "-b:a:c":
      ibrav_ = -9 
    elif alt_axes == "bcoA-type":
      ibrav_ = 91 
  elif bravais_index == 12 or bravais_index == 13:
    if alt_axes == "unique-axis-b":
      ibrav_ = -bravais_index 
  else:
    ibrav_ = bravais_index
  return ibrav_ 

def at2celldm(at, ibrav_, alat_ = None):
  c1 = np.sqrt(at[0].dot(at[0]))
  c2 = None
  c3 = None
  c4 = None
  c5 = None
  c6 = None
  if ibrav_ == 0:
    if alat_:
      c1 = alat_ 
  elif ibrav_  == 2:
    c1 = c1 * np.sqrt(2.e0)
  elif ibrav_ == 3 or ibrav_ == -3:
    c1 = c1 / np.sqrt(3.e0)*2.e0  
  elif ibrav_ == 4 or ibrav_ == 6:
    c3 = np.sqrt(at[2].dot(at[2]))/c1
  elif ibrav_ == 5 or ibrav_ == -5:
    c4 = at[0].dot(at[1])/c1/np.sqrt(at[1].dot(at[1]))
  elif ibrav_ == 7:
    c1 = np.abs(at[0][0]) * 2.e0 
    c3 = np.abs(at[0][2]/at[0][0])
  elif ibrav_ == 8:
    c2 = np.sqrt(at[1].dot(at[1]))/c1 
    c3 = np.sqrt(at[2].dot(at[2]))/c1
  elif ibrav_ == 9 or ibrav_ == -9:
    c1 = np.abs(at[0][0]) * 2.e0
    c2 = np.abs(at[1][1]) * 2.e0 / c1 
    c3 = np.abs(at[2][2]) / c1 
  elif ibrav_ == 91:
    c2 = np.abs(at[1][1]) * 2.e0 / c1 
    c3 = np.abs(at[2][2]) * 2.e0 / c1 
  elif ibrav_ == 10:
    c1 = np.abs(at[0][0]) * 2.e0
    c2 = np.abs(at[1][1]) * 2.e0 / c1 
    c3 = np.abs(at[2][2]) * 2.e0 / c1
  elif ibrav_ == 11:
    c1 = np.abs(at[0][0]) * 2.e0
    c2 = np.abs(at[0][1]) * 2.e0 / c1 
    c3 = np.abs(at[0][2]) * 2.e0 / c1 
  elif ibrav_ == 12 or ibrav_ == -12:
    c2 = np.sqrt(at[1].dot(at[1])) / c1 
    c3 = np.sqrt(at[2].dot(at[2])) / c1
    if ibrav_  == 12:
      c4 = at[0].dot(at[1]) /c1 / c1 / c2 
    else: 
      c5 =  at[0].dot(at[2]) /c1 / c1 / c3
  elif ibrav_ == 13:
     c1 = np.abs(at[0][0]) * 2.e0 
     c2 = np.sqrt(at[1].dot(at[1])) / c1 
     c3 = np.abs(at[0][2]/at[0][0])
     c4 = at[1][0]/at[0][0] /c2 / 2.e0
  elif ibrav_ == -13:
    c1 = np.abs(at[0][0]) *2.e0 
    c2 = np.abs (at[1][1]/at[1][0])
    c3 = np.sqrt(at[2].dot(at[3]))/ c1 
    c5 = at[2][0]/ at[0][0]/ c3/ 2.e0
  elif ibrav_ == 14:
    c2 = np.sqrt(at[1].dot(at[1]))/c1 
    c3 = np.sqrt(at[2].dot(at[2]))/c1
    c4 = at[1].dot(at[2])/ c1 / c1 / c2 / c3
    c5 = at[0].dot(at[2])/ c1 / c1 / c3
    c6 = at[1].dot(at[0])/ c1 / c1 / c2 
  return (c1,c2,c3,c4,c5,c6)

def qe_abc2at(abc, ibrav_, abc_is_primitive = False):
  celldm_ =  abc2celldm(abc, ibrav_, abc_is_primitive)
  at_ = latgen(celldm, ibrav_)





class lattice:
  def __init__(self, vectors = None, abc = None, abc_is_primitive = False, 
               bravais_index = None, alt_axes = None, celldm = None, alat_ = None):
    self.bravais_index = bravais_index
    self.alt_axes =of alt_axes 
    self.cell  = self.__init_vectors(vectors, abc, abc_is_primitive, celldm)
    self.ibrav, self.celldm = self.get_ibrav_celldm() 
    abc_ = self.get_primitive_abc(vectors, abc, abc_is_primitive)
    self.primitive_abc = dict (zip(("a","b","c","cosBC","cosCA","cosAB"),abc_)) 
    abc_ = self.get_conventional_abc()
    self.conventional_abc = dict(zip(("a","b","c","cosBC","cosCA","cosAB")), abc_)
    self.alt_axes = alt_axes
    
  def __init_vectors(self, vectors = None, abc = None, abc_is_primitive = False, celldm = None):
    ibrav_ = qe_ibrav(self.bravais_index, self.alt_axes)
    if vectors:
      cellvec = np.asarray(vectors)
    elif abc:
      cellvec = qe_abc2at(abc, abc_is_primitive, ibrav_)
    elif celldm:
      cellvec = latgen(ibrav_, celldm)
    else:
      logger.error("lattice initialization failed")
    return cellvec 

  def get_ibrav_celldm(self, alat_ = None):
    ibrav_ = qe_ibrav(self.bravais_index, self.alt_axes)
    celldm = at2celldm(self.cell, ibrav_) 
    return ibrav_, celldm 
  
  def get_primitive_abc(abc = None, abc_is_primitive = False):
    if abc:
      if (np.round(self.cell - qe_abc2at(abc, self.ibrav, abc_is_primitive), 5) != 0.0).any(): 
          logger.error("abc parameters in input do not match with expected cell vectors")  
    a = round(np.sqrt(self.cell[0].dot(self.cell[0])),5) 
    b = round(np.sqrt(self.cell[1].dot(self.cell[1])),5)
    c = round(np.sqrt(self.cell[1].dot(self.cell[1])),5)
    bc = round(self.cell[1].dot(self.cell[2])/b/c, 5)   
    ca = round(self.cell[2].dot(self.cell[0])/c/a, 5)
    ab = round(self.cell[0].dot(self.cell[1])/a/b, 5)
    return a, b, c, bc, ca, ab

  def get_conventional_abc(self):
    abc = qe_at2abc(self.cell, self.ibrav)
    return abc 











