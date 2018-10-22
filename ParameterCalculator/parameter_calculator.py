import numpy as np
import logging
from astropy import constants as const

M_SUN = const.M_sun
G = const.G
C = const.c

# Sample default parameter values for M,S,L
_M1 = 1.5 * M_SUN
_M2 = 1.5 * M_SUN
_S1 = np.array([2e42, 0., 0.])
_S2 = np.array([3e42, 10.5, 0])
_L = np.array([4e43, 10.5, 10.5])

class ParameterCalculator:
    def __init__(self):
        self.m1 = _M1
        self.m2 = _M2
        self.S1 = _S1
        self.S2 = _S2
        self.L = _L

        self.logger = logging.getLogger(__name__)
        self.logger.info('Created ParameterCalculator with default parameters\n{}'.format(self))

    def __str__(self):
        return 'Mass1: {}, Mass2: {}, S1: {}, S2: {}, L: {}'.format(self.m1, self.m2, self.S1, self.S2, self.L)

    @property
    def M(self):
        return self.m1 + self.m2

    @property
    def b1(self):
        return 2 + 3 * self.q() / 2

    @property
    def b2(self):
        return 2 + 2 * self.q() / 3

    def q(self):
        '''Return the Mass Ratio'''
        return self.m1 / self.m2

    def chirpMass(self):
        '''Return the Chirp Mass'''
        return (self.m1 * self.m2)**(3./5) / self.M**(1./5)

    def eta(self):
        '''Return the Symmetric Mass Ratio eta'''
        return (self.m1 * self.m2) / self.M

    def chi_eff(self):
        '''
        Return the Effective Spin Parameter chi_eff
        (component spins in the aligned direction)
        '''
        return ((C / G / self.M) * 
                np.dot((self.S1 / self.m1 + self.S2 / self.m2),
                       (self.L / np.linalg.norm(self.L))))

    def chi_p(self):
        '''
        Return the Precession Spin Parameter chi_p
        (component in unaligned direction)
        '''
        s1_perp = np.linalg.norm(np.cross(self.S1, self.L) / np.linalg.norm(self.L))
        s2_perp = np.linalg.norm(np.cross(self.S2, self.L) / np.linalg.norm(self.L))
        return ((C / (self.b1 * G * self.m1**2)) * 
                max(self.b1 * s1_perp, self.b2 * s2_perp))
