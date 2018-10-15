#!/usr/local/bin/python3

from ParameterCalculator.parameter_calculator import ParameterCalculator, SOLAR_MASS

import numpy as np

calc = ParameterCalculator()
calc.m1 = 3.5 * SOLAR_MASS

print('q:\t\t', calc.q())
print('Chirp Mass:\t', calc.chirpMass())
print('eta:\t\t', calc.eta())
print('Chi effective:\t', calc.chi_eff())
print('Chi P:\t\t', calc.chi_p())
