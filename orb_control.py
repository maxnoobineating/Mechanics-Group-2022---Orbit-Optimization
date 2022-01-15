import numpy as np
from math import sqrt, cos, sin, acos, exp, log, pi
from itertools import islice
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy import optimize as opt
import scipy.integrate as itg
import scipy.interpolate as itp
from matplotlib import pyplot as plt
from numba import njit, float64, int64
import cProfile
import pstats
import copy

# All the constant, state variable, and imports are all here

# getter always return with scaled factor add-ons (decorated constant)
# , except itself!

# Respect My Constraint!!!!!!!!!!!!!!
# for weighting the constraint & objectives
Ratio = 1e3

class Control:
    # Scaling factors
    # choose the scale level so that the number
    # fed into optimizer is closer to 1
    _KM = 1e6 # 1e6 km
    _S = 24*3600 # 1 day of seconds

    @property
    def KM(self):
        self.broadcast()
        return self._KM
    @property
    def S(self):
        self.broadcast()
        return self._S

    # Switch on/off the scale factors
    _switch = False

    def toScale(self, x):
        if self._switch:
            return 1
        else:
            return x

    def switch(self, switch_to=None):
        if switch == None:
            self._switch = not(self._switch)
        else:
            return switch_to

    # Boadcast & dumper
    def broadcast(self):
        pass
        # print()
        # print("Oops! Be sure that we're not running an optimization!!")
#        print("args | mu: {0:.3f}| aE: {1:.3f}| eE: {2:.3f}| aM: {3:.3f}|\
# eM: {4:.3f}| periM: {5:.3f}| Isp: {6:.3f}| fM: {7:.3f}| DV: {8:.3f}|"\
#            .format(self.mu, self.aE, self.eE, self.aM, self.eM, self.periM, self.Isp, self.fM, self.DV))

    # def dumpAll():


    # Minimizer parameter
    @property
    def xtol(self):
        return self._xtol
    @xtol.setter
    def xtol(self, a):
        self._xtol = a

    @property
    def gtol(self):
        return self._gtol
    @gtol.setter
    def gtol(self, a):
        self._gtol = a

    @property
    def maxiter(self):
        return self._maxiter
    @maxiter.setter
    def maxiter(self, a):
        self._maxiter = a

    @property
    def verb(self):
        return self._gtol
    @verb.setter
    def verb(self, a):
        self._gtol = a


    # Iteration Parameter
    # fid = 2**nfld+1
    @property
    def nfld(self):
        return self._nfld
    @nfld.setter
    def nfld(self, a):
        print(a, "folds, Fidelity set: ", (2**a)+1)
        self._nfld = a

    @property
    def fid(self):
        return 2**self._nfld+1


    # Orbit Elements & Spacecraft specs
    # KM^3/S^2
    @property
    def mu(self):
        self.broadcast()
        return self._mu/self.toScale(self._KM**3/self._S**2)
    @mu.setter
    def mu(self, a):
        self._mu = a
        self.broadcast()

    # KM
    @property
    def aE(self):
        self.broadcast()
        return self._aE/self.toScale(self._KM)
    @aE.setter
    def aE(self, a):
        self._aE = a
        self.broadcast()

    @property
    def eE(self):
        self.broadcast()
        return self._eE
    @eE.setter
    def eE(self, a):
        self._eE = a
        self.broadcast()

    @property
    def periE(self):
        self.broadcast()
        return self._periE
    @periE.setter
    def periE(self, a):
        self._periE = a
        self.broadcast()

    # KM
    @property
    def aM(self):
        self.broadcast()
        return self._aM/self.toScale(self._KM)
    @aM.setter
    def aM(self, a):
        self._aM = a
        self.broadcast()

    @property
    def eM(self):
        self.broadcast()
        return self._eM
    @eM.setter
    def eM(self, a):
        self._eM = a
        self.broadcast()

    @property
    def periM(self):
        self.broadcast()
        return self._periM
    @periM.setter
    def periM(self, a):
        self._periM = a
        self.broadcast()

    @property
    def InitMass(self):
        self.broadcast()
        return self._InitMass
    @InitMass.setter
    def InitMass(self, a):
        self._InitMass = a
        self.broadcast()

    # KM/S
    @property
    def Isp(self):
        self.broadcast()
        return self._Isp*9.80665/1000/self.toScale(self._KM/self._S)
    @Isp.setter
    def Isp(self, a):
        self._Isp = a
        self.broadcast()

    # KM/S^2
    @property
    def fM(self):
        self.broadcast()
        return self._fM/self.InitMass/self.toScale(self._KM/self._S**2)
    @fM.setter
    def fM(self, a):
        self._fM = a
        self.broadcast()

    @property
    def Mr(self):
        self.broadcast()
        return self._Mr
    @Mr.setter
    def Mr(self, a):
        self._Mr = a
        self.broadcast()

    # KM/S
    @property
    def DV(self):
        self.broadcast()
        return self.Isp*log(1/self.Mr)


    # Generated contents
    # Control Orbit (plot & optimization parameter)
    @property
    def elpsE(self):
        self.broadcast()
        return {'a':self.aE, 'e':self.eE, 'peri':self.periE, 'fs': (0, 2*pi)}

    @property
    def elpsM(self):
        self.broadcast()
        return {'a':self.aM, 'e':self.eM, 'peri':self.periM, 'fs': (0, 2*pi)}

    # Boundary Condition & constraint
    @property
    def rslb(self):
        self.broadcast()
        return [0 for i in range(self.fid)]
    @property
    def cslb(self):
        self.broadcast()
        return [-np.inf for i in range(self.fid-1)]
    @property
    def thslb(self):
        self.broadcast()
        return [-np.inf, -np.inf]
    @property
    def lbs(self):
        self.broadcast()
        return self.rslb+self.cslb+self.thslb

    @property
    def rsub(self):
        self.broadcast()
        return [np.inf for i in range(self.fid)]
    @property
    def csub(self):
        self.broadcast()
        return [np.inf for i in range(self.fid-1)]
    @property
    def thsub(self):
        self.broadcast()
        return [np.inf, np.inf]
    @property
    def ubs(self):
        self.broadcast()
        return self.rsub+self.csub+self.thsub

    @property
    def args(self):
        self.broadcast()
        return (self.mu, self.aE, self.eE, self.aM, self.eM, self.periM, self.Isp, self.fM, self.DV)


#--------------------------------------------#
Earth_Mars = Control()

Earth_Mars.xtol = 1e-2
Earth_Mars.gtol = 1e-8
Earth_Mars.maxiter = 5e2
Earth_Mars.verb = 0

Earth_Mars.nfld = 5
Earth_Mars.fid # == 1*(2**nfld)+1

Earth_Mars.mu = 132712*1e6 # km^3/s^2

# Earth
# absolute anamoly set to be Earth's anamoly
Earth_Mars.aE = 149.598*1e6 # km
Earth_Mars.eE = 0.0167
Earth_Mars.periE = 0.0*pi

# Mars
Earth_Mars.aM = 227.956*1e6 # km
Earth_Mars.eM = 0.0935
Earth_Mars.periM = (336.04084-102.94719)*pi/360 # rad

# Control Orbit (plot & optimization parameter)
Earth_Mars.elpsE # == {'a':aE, 'e':eE, 'peri':periE, 'fs': (0, 2*pi)}
Earth_Mars.elpsM # == {'a':aM, 'e':eM, 'peri':periM, 'fs': (0, 2*pi)}

# Boundary Condition & constraint
Earth_Mars.rslb # == [1e-5 for i in range(fid)]
Earth_Mars.cslb # == [-np.inf for i in range(fid-1)]
Earth_Mars.thslb # == [-np.inf, -np.inf]
Earth_Mars.lbs # == rslb+cslb+thslb

Earth_Mars.rsub # == [np.inf for i in range(fid)]
Earth_Mars.csub # == [np.inf for i in range(fid-1)]
Earth_Mars.thsub # == [np.inf, np.inf]
Earth_Mars.ubs # == rsub+csub+thsub

EM = {}
#--------------------------------------------#
EM['RD-0410'] = copy.deepcopy(Earth_Mars)

# Spacecraft
# 1. EM['RD-0410']
# RD-0410 thermal nuclear
# Isp = 910 s
# thrust = 35.30 kN/vheicle initial mass (1e3*m*kg/s^2)
# Isp = 9.80665*910/1000 # = 9.80665*x*1e-3 (km/s), x be the specific impulse measured by weight
# fM = 35.3/5000 # km*kg/s^2

EM['RD-0410'].InitMass = 50e3 # kg
EM['RD-0410'].Isp = 5000 # = 9.80665*x*1e-3 (km/s), x be the specific impulse measured by weight
EM['RD-0410'].fM = 35.5e-3 # (km*kg/s^2)
EM['RD-0410'].Mr = 0.1 # dry/full, for constraint, Option B
EM['RD-0410'].DV # == Isp*log(1/Mr) # km/s

EM['RD-0410'].args # == (mu, aE, eE, aM, eM, periM, Isp, fM, DV)

#--------------------------------------------#

EM['VASIMR'] = copy.deepcopy(Earth_Mars)

# Spacecraft
# 2. EM['VASIMR']
# VASIMR
# Isp = 5000
# thrust = 5 kN

EM['VASIMR'].InitMass = 50e3 # kg
EM['VASIMR'].Isp = 5000 # = 9.80665*x*1e-3 (km/s), x be the specific impulse measured by weight
EM['VASIMR'].fM = 5e-3 # (km*kg/s^2)
EM['VASIMR'].Mr = 0.1 # dry/full, for constraint, Option B
EM['VASIMR'].DV # == Isp*log(1/Mr) # km/s

EM['VASIMR'].args # == (mu, aE, eE, aM, eM, periM, Isp, fM, DV)

#--------------------------------------------#
