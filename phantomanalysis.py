#
# Analyse a Phantom dump file with python
#
# Written by David Liptai 2018
#

import pyphantom
from ctypes import *
import numpy as np
import matplotlib.pyplot as plt

class PhantomAnalysis(pyphantom.Simulation):
    def __init__(self, dumpfile):
      self.libph = CDLL('./libphantom.so')

      # Read dumpfile
      time,hfact,massofgas,massofdust = self.read_dump(dumpfile)

      # Units
      udist, umass, utime, udens, umagfd = self.get_units()
      self.units = {'udist':udist,
                    'umass':umass,
                    'utime':utime,
                    'udens':udens,
                    'umagfd':umagfd}

      # Note that we divide by code units so that quantities are in code units by default

      # General stuff
      self.time   = time
      self.hfact  = hfact
      self.massofgas  = massofgas/umass
      self.massofdust = massofdust/umass

      # Gas quantities
      npart = self.get_npart()
      self.npart       = npart
      self.xyzh        = self.get_part_xyzh(npart)/udist
  
      #self.vxyz        = self.get_part_vxyz(npart)/(udist/utime)

      # Removed by cadman@roe.ac.uk as these values are not stored in growingdisc simulations
      #self.utherm      = self.get_part_u(npart)/(udist**2/utime**2)
      #self.temperature = self.get_part_temp(npart)
      #self.bxyz        = self.get_part_bxyz(npart)/umagfd

      self.grainsize   = self.get_part_grainsize(npart)

      # Point masses
      nptmass             = self.get_nptmass()
      self.nptmass        = nptmass
      self.ptmass_xyzmh   = self.get_ptmass_xyzmh(nptmass)
      self.ptmass_xyzmh[0,:] = self.ptmass_xyzmh[0,:]/udist
      self.ptmass_xyzmh[1,:] = self.ptmass_xyzmh[1,:]/udist
      self.ptmass_xyzmh[2,:] = self.ptmass_xyzmh[2,:]/udist
      self.ptmass_xyzmh[3,:] = self.ptmass_xyzmh[3,:]/umass
      self.ptmass_xyzmh[4,:] = self.ptmass_xyzmh[4,:]/udist
      self.ptmass_vxyz    = self.get_ptmass_vxyz(nptmass)/(udist/utime)
      self.ptmass_spinxyz = self.get_ptmass_spinxyz(nptmass)/(udist**2*umass/utime)

      self.itype, self.dustfrac, self.utherm = self.get_extravals(npart)
      self.dustfrac
      self.density = np.zeros(npart)
      self.density[self.itype[:]==1] = self.massofgas/(abs(self.xyzh[3][self.itype==1])**3)
      self.density[self.itype[:]==2] = self.massofdust/(abs(self.xyzh[3][self.itype==2])**3)


