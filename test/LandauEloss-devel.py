#!/usr/bin/env python
# coding: utf-8

#%pylab
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import matplotlib.dates as dts
import numpy as np
import pandas as pd
import itertools
import os
import datetime
from scipy.optimize import curve_fit
from scipy.misc import factorial
plt.rcParams.update({'font.size': 18})

from ELOSS import *

KEvals = np.linspace(10.,100000.,10000)
dpvals04 = []
dpvals03 = []
dpvals05 = []
dpvals30 = []
dedxvals = []

dpsmall = []

x = 0.49 # cm

for KE in KEvals:
    dpvals04.append(dpdx(KE,0.4,Mmu))
    dpvals05.append(dpdx(KE,0.5,Mmu))
    dpvals03.append(dpdx(KE,0.3,Mmu))
    dpvals30.append(dpdx(KE,3.0,Mmu))
    dpsmall.append(dpdx(KE,0.03,Mmu))
    dedxvals.append(rho * dedx(KE,Mmu,dens=True))

fig = plt.figure(figsize=(6,6))
plt.plot(KEvals,dpvals05,'-',color='b' ,lw=3,label='$\Delta p/dx$ @ 0.5 cm')
plt.plot(KEvals,dpvals03,'--',color='b',lw=3,label='$\Delta p/dx$ @ 0.3 cm')
plt.plot(KEvals,dpvals30,'-.',color='b',lw=3,label='$\Delta p/dx$ @ 3.0 cm')
plt.plot(KEvals,dedxvals,color='r',lw=2,label='$<dE/dx>$')
plt.grid(which='both')
plt.xlabel('Muon Kinetic Energy [MeV]',fontsize=16,fontweight='bold')
plt.ylabel('MeV/cm',fontsize=16,fontweight='bold')
#plt.ylim([1.75,1.85])
plt.ylim([1.45,3.12])
plt.legend(loc=1,fontsize=16)
plt.xlim([50,100000])
plt.xscale('log')
plt.title('Collision Stopping Power for $\mu$ in LAr',fontsize=16,fontweight='bold')
plt.show()

