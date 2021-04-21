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


def abc_from_cell(cell, dgts=5):
  """
  cell: dictionary with 'a1','a2','a3' keys with values corresponding to the 3 lattice vectors
  returns: 6ple with values of A,B,C,cos(AB), cos(BC), cos(AC) lenghts in Angstrom and cosines of the 
           lattice vectors
  """
  at = np.array((cell['a1'],cell['a2'], cell['a3']))
  A = round(np.sqrt(at[0].dot(at[0])),dgts) 
  B = round(np.sqrt(at[1].dot(at[1])),dgts)
  C = round(np.sqrt(at[2].dot(at[2])),dgts) 
  COSAB = round(at[0].dot(at[1])/A/B, dgts) 
  COSBC = round(at[1].dot(at[2])/B/C, dgts)
  COSAC = round(at[2].dot(at[0])/C/A, dgts)  
  bohr2ang = 0.529177 
  return round(A * bohr0ang,dgts), round(B * bohr2ang), round(C * bohr2ang,dgts),
     COSAB, COSBC, COSAC 



  bi_ = ibrav_from_abc(A, B, C, COSAB, COSBC, COSAC)
  if isinstance(bi_, tuple):
    if bi_ == (5,-5):
      bi = resolve_trigonal(at)
    elif bi_ == (7,3):
      bi = resolve_tbc(at)
    elif bi_ == (9,-9,14):
      bi = resolve_bco(at)
    elif bi_ == (11, 14):
      bi = resolve_obc
    elif bi_ == (10,14):
      bi = resolve_of(at)
  else:
    bi = bi_ 
  

def ibrav_from_abc(a,b,c,cab,cbc,cac,dgts=6):
  """
  a,b,c,cab,cbc,cac:  lenghts of lattice vectors and cosines of the angles they span.
  dgts: integer chose precision for numerical checks.   
  """
  def eqq(f,g,dgts):
    return isclose (f,g,abs_tol=0.5*10^(-dgts))
  def eqq0(f,dgts):
    return isclose(f,0.,rel_tol=0.5*10^(-dgts))
  if eqq(a,b) and eqq(b,c):
    if eqq(cab,cac) and eqq(cab,cbcb):
      if eqq0(cab):
        return 1
      elif eqq(cab,0.5):
        return 2
      elif eqq(cab, -1.0/3.0):
        return -3 
      else:
        return (5,-5)
    elif eqq (cab, cac) and not eqq(cab,cbc):
      return (7,3)
    elif eqq(cab,-cac) and eqq(cab, cbc) and eqq(cab, 1./3.): 
      return 3
    else:
      return (11,14) 
  elif eqq(a,b) and not eqq(a,c):
    if eqq0(cab) and eqq0(cac) and eqq0(cbc):
      return 6
    elif eqq(cab,-0.5) and eqq0(cac) and eqq0(cbc):
      return 4
    elif eqq0(cac) and eqq0(cbc):
      return (9,-9,14)  
    elif eqq(cac, -cbc):
      return -13 
  elif eqq(a,c) and not eqq(a,b):
    return 13 
  elif eqq(b,c) and not eqq(a,b):
    return 91
  elif not eqq(a,b) and not eqq(a,c) and not eqq(b,c):
    if eqq0(cab) and eqq0(cac) and eqq0(cbc):
      return 8
    elif not eqq0(cab) and eqq0(cac) and eqq0(cbc):
      return 12
    elif eqq0(cab) and not eqq0(cac) and eqq0(cbc):
      return -12
    elif not eqq0(cab) and not eqq0(cac) and not eqq0(cbc):
      return (10,14) 
  return 14

def resolve_trigonal(at,dgts=6):
  def eqq(f,g):
    return isclose(f, g, abs_tol=0.5*10^(-dgts))
  if eqq(at[0,2],at[1,2]) and eqq(at[1,2],at[2,2]):
    return 5
  else:
    return -5

def resolve_tbc(at,dgts=6):
  def eqq(f,g):
    return isclose(f,g,abs_tol=0.5*10^(-dgts))
  if eqq(at[0,0],at[0,1]) and eqq(at[1,0],at[1,1]):
    return 7
  else:
    return 3

def resolve_bco(at, dgts=6):
  def eqq(f,g):
    return isclose(f,g,abs_tol=0.5*10^(-dgts))
  if eqq(at[0,0],at[1,0]) and eqq(at[0,1],-at[1,1]):
    return -9 
  elif eqq(at[0,0],-at[1,0]) and eqq(at[0,1],at[1,1]):
    return 9
  else:
    return 14

def resolve_obc(at, dgts=6):
  def eqqabs(f,g):
    return isclose (abs(f),abs(g),abs_tol=0.5*10^(-dgts))
  if eqqabs(at[0,0], at[1,0]) and eqqabs(at[0,1],at[1,1]):
    return 11 
  else:
    return 14 
  

def resolve_of(at, dgts=6):
  def eqq(f,g):
    return isclose(f,g,abs_tol=0.5*10^(-dgts))
  if 