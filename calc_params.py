#!/usr/local/bin/python3

from ParameterCalculator.parameter_calculator import ParameterCalculator, M_SUN

import numpy as np


calc = ParameterCalculator()
calc.m1 = 3.5 * M_SUN

print('q:\t\t', calc.q())
print('Chirp Mass:\t', calc.chirpMass())
print('eta:\t\t', calc.eta())
print('Chi effective:\t', calc.chi_eff())
print('Chi P:\t\t', calc.chi_p())
