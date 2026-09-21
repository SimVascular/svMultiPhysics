# SPDX-License-Identifier: CC-BY-3.0
# Adapted from Physiome CellML API-generated Python for FitzHugh (1961),
# "Impulses and Physiological States in Theoretical Models of Nerve Membrane".
# CellML model author: Penny Noble; source documentation: Catherine Lloyd.
# License: https://creativecommons.org/licenses/by/3.0/
# Changes: caller-supplied stimulus; unused generated helpers, imports,
# demonstration solver, and plotting removed. Rate expressions are unchanged.
# See README.md for the exact source, revision, upstream hash, and attribution.

# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 2
sizeConstants = 3

# States: v, w. Constants: alpha, gamma, epsilon. Time is in milliseconds.
def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0
    states[1] = 0
    constants[0] = -0.08
    constants[1] = 3
    constants[2] = 0.005
    return (states, constants)

def computeRates(voi, states, constants, stimulus):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = 1.00000*constants[2]*(states[0]-constants[1]*states[1])
    algebraic[0] = stimulus
    rates[0] = 1.00000*((states[0]*(states[0]-constants[0])*(1.00000-states[0])-states[1])+algebraic[0])
    return(rates)
