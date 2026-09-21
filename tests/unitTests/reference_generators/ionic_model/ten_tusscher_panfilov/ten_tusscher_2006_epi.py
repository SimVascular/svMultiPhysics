# SPDX-License-Identifier: CC-BY-3.0
# Adapted from Physiome CellML API-generated Python for the TP06 epicardial
# model (ten_tusscher_model_2006_IK1Ko_epi_units.cellml).
# CellML model author: Penny Noble.
# License: https://creativecommons.org/licenses/by/3.0/
# Changes: caller-supplied stimulus; algebraic workspace returned for
# Rush--Larsen updates; unused generated legends, vectorized output, imports,
# demonstration solver, and plotting removed. NumPy scalar helpers were
# replaced with equivalent standard-library operations. Scientific constants
# and all other generated rate expressions are unchanged.
# See README.md for the exact source, revision, upstream hash, and attribution.

# Size of variable arrays:
sizeAlgebraic = 70
sizeStates = 19
sizeConstants = 53

from math import exp, log


def power(base, exponent):
    return base ** exponent


def less(left, right):
    return left < right


def custom_piecewise(cases):
    for condition, value in zip(cases[0::2], cases[1::2]):
        if condition:
            return value
    raise ValueError("piecewise expression has no matching case")


def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -85.23
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 0.185
    constants[4] = 0.016404
    constants[5] = 10
    constants[6] = 1000
    constants[7] = 1
    constants[8] = 52
    constants[9] = 0.03
    constants[10] = 5.4
    constants[11] = 140
    states[1] = 136.89
    states[2] = 8.604
    constants[12] = 2
    states[3] = 0.000126
    constants[13] = 5.405
    constants[14] = 0.153
    states[4] = 0.00621
    states[5] = 0.4712
    constants[15] = 0.392
    states[6] = 0.0095
    constants[16] = 14.838
    states[7] = 0.00172
    states[8] = 0.7444
    states[9] = 0.7045
    constants[17] = 0.00029
    constants[18] = 0.0000398
    states[10] = 0.00036
    states[11] = 3.373e-5
    states[12] = 0.7888
    states[13] = 0.9755
    states[14] = 0.9953
    constants[19] = 0.000592
    constants[20] = 0.294
    states[15] = 0.999998
    states[16] = 2.42e-8
    constants[21] = 2.724
    constants[22] = 1
    constants[23] = 40
    constants[24] = 1000
    constants[25] = 0.1
    constants[26] = 2.5
    constants[27] = 0.35
    constants[28] = 1.38
    constants[29] = 87.5
    constants[30] = 0.1238
    constants[31] = 0.0005
    constants[32] = 0.0146
    states[17] = 3.64
    states[18] = 0.9073
    constants[33] = 0.15
    constants[34] = 0.045
    constants[35] = 0.06
    constants[36] = 0.005
    constants[37] = 1.5
    constants[38] = 2.5
    constants[39] = 1
    constants[40] = 0.102
    constants[41] = 0.0038
    constants[42] = 0.00025
    constants[43] = 0.00036
    constants[44] = 0.006375
    constants[45] = 0.2
    constants[46] = 0.001
    constants[47] = 10
    constants[48] = 0.3
    constants[49] = 0.4
    constants[50] = 0.00025
    constants[51] = 0.001094
    constants[52] = 0.00005468
    return (states, constants)

def computeRates(voi, states, constants, stimulus):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+20.0000)/7.00000))
    algebraic[20] = 1102.50*exp(-(power(states[0]+27.0000, 2.00000))/225.000)+200.000/(1.00000+exp((13.0000-states[0])/10.0000))+180.000/(1.00000+exp((states[0]+30.0000)/10.0000))+20.0000
    rates[12] = (algebraic[7]-states[12])/algebraic[20]
    algebraic[8] = 0.670000/(1.00000+exp((states[0]+35.0000)/7.00000))+0.330000
    algebraic[21] = 562.000*exp(-(power(states[0]+27.0000, 2.00000))/240.000)+31.0000/(1.00000+exp((25.0000-states[0])/10.0000))+80.0000/(1.00000+exp((states[0]+30.0000)/10.0000))
    rates[13] = (algebraic[8]-states[13])/algebraic[21]
    algebraic[9] = 0.600000/(1.00000+power(states[10]/0.0500000, 2.00000))+0.400000
    algebraic[22] = 80.0000/(1.00000+power(states[10]/0.0500000, 2.00000))+2.00000
    rates[14] = (algebraic[9]-states[14])/algebraic[22]
    algebraic[10] = 1.00000/(1.00000+exp((states[0]+20.0000)/5.00000))
    algebraic[23] = 85.0000*exp(-(power(states[0]+45.0000, 2.00000))/320.000)+5.00000/(1.00000+exp((states[0]-20.0000)/5.00000))+3.00000
    rates[15] = (algebraic[10]-states[15])/algebraic[23]
    algebraic[11] = 1.00000/(1.00000+exp((20.0000-states[0])/6.00000))
    algebraic[24] = 9.50000*exp(-(power(states[0]+40.0000, 2.00000))/1800.00)+0.800000
    rates[16] = (algebraic[11]-states[16])/algebraic[24]
    algebraic[0] = 1.00000/(1.00000+exp((-26.0000-states[0])/7.00000))
    algebraic[13] = 450.000/(1.00000+exp((-45.0000-states[0])/10.0000))
    algebraic[26] = 6.00000/(1.00000+exp((states[0]+30.0000)/11.5000))
    algebraic[34] = 1.00000*algebraic[13]*algebraic[26]
    rates[4] = (algebraic[0]-states[4])/algebraic[34]
    algebraic[1] = 1.00000/(1.00000+exp((states[0]+88.0000)/24.0000))
    algebraic[14] = 3.00000/(1.00000+exp((-60.0000-states[0])/20.0000))
    algebraic[27] = 1.12000/(1.00000+exp((states[0]-60.0000)/20.0000))
    algebraic[35] = 1.00000*algebraic[14]*algebraic[27]
    rates[5] = (algebraic[1]-states[5])/algebraic[35]
    algebraic[2] = 1.00000/(1.00000+exp((-5.00000-states[0])/14.0000))
    algebraic[15] = 1400.00/(power(1.00000+exp((5.00000-states[0])/6.00000), 1.0/2))
    algebraic[28] = 1.00000/(1.00000+exp((states[0]-35.0000)/15.0000))
    algebraic[36] = 1.00000*algebraic[15]*algebraic[28]+80.0000
    rates[6] = (algebraic[2]-states[6])/algebraic[36]
    algebraic[3] = 1.00000/(power(1.00000+exp((-56.8600-states[0])/9.03000), 2.00000))
    algebraic[16] = 1.00000/(1.00000+exp((-60.0000-states[0])/5.00000))
    algebraic[29] = 0.100000/(1.00000+exp((states[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((states[0]-50.0000)/200.000))
    algebraic[37] = 1.00000*algebraic[16]*algebraic[29]
    rates[7] = (algebraic[3]-states[7])/algebraic[37]
    algebraic[4] = 1.00000/(power(1.00000+exp((states[0]+71.5500)/7.43000), 2.00000))
    algebraic[17] = custom_piecewise([less(states[0] , -40.0000), 0.0570000*exp(-(states[0]+80.0000)/6.80000) , True, 0.00000])
    algebraic[30] = custom_piecewise([less(states[0] , -40.0000), 2.70000*exp(0.0790000*states[0])+310000.*exp(0.348500*states[0]) , True, 0.770000/(0.130000*(1.00000+exp((states[0]+10.6600)/-11.1000)))])
    algebraic[38] = 1.00000/(algebraic[17]+algebraic[30])
    rates[8] = (algebraic[4]-states[8])/algebraic[38]
    algebraic[5] = 1.00000/(power(1.00000+exp((states[0]+71.5500)/7.43000), 2.00000))
    algebraic[18] = custom_piecewise([less(states[0] , -40.0000), (((-25428.0*exp(0.244400*states[0])-6.94800e-06*exp(-0.0439100*states[0]))*(states[0]+37.7800))/1.00000)/(1.00000+exp(0.311000*(states[0]+79.2300))) , True, 0.00000])
    algebraic[31] = custom_piecewise([less(states[0] , -40.0000), (0.0242400*exp(-0.0105200*states[0]))/(1.00000+exp(-0.137800*(states[0]+40.1400))) , True, (0.600000*exp(0.0570000*states[0]))/(1.00000+exp(-0.100000*(states[0]+32.0000)))])
    algebraic[39] = 1.00000/(algebraic[18]+algebraic[31])
    rates[9] = (algebraic[5]-states[9])/algebraic[39]
    algebraic[6] = 1.00000/(1.00000+exp((-8.00000-states[0])/7.50000))
    algebraic[19] = 1.40000/(1.00000+exp((-35.0000-states[0])/13.0000))+0.250000
    algebraic[32] = 1.40000/(1.00000+exp((states[0]+5.00000)/5.00000))
    algebraic[40] = 1.00000/(1.00000+exp((50.0000-states[0])/20.0000))
    algebraic[42] = 1.00000*algebraic[19]*algebraic[32]+algebraic[40]
    rates[11] = (algebraic[6]-states[11])/algebraic[42]
    algebraic[55] = ((((constants[21]*constants[10])/(constants[10]+constants[22]))*states[2])/(states[2]+constants[23]))/(1.00000+0.124500*exp((-0.100000*states[0]*constants[2])/(constants[0]*constants[1]))+0.0353000*exp((-states[0]*constants[2])/(constants[0]*constants[1])))
    algebraic[25] = ((constants[0]*constants[1])/constants[2])*log(constants[11]/states[2])
    algebraic[50] = constants[16]*(power(states[7], 3.00000))*states[8]*states[9]*(states[0]-algebraic[25])
    algebraic[51] = constants[17]*(states[0]-algebraic[25])
    algebraic[56] = (constants[24]*(exp((constants[27]*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], 3.00000))*constants[12]-exp(((constants[27]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[11], 3.00000))*states[3]*constants[26]))/((power(constants[29], 3.00000)+power(constants[11], 3.00000))*(constants[28]+constants[12])*(1.00000+constants[25]*exp(((constants[27]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))))
    rates[2] = ((-1.00000*(algebraic[50]+algebraic[51]+3.00000*algebraic[55]+3.00000*algebraic[56]))/(1.00000*constants[4]*constants[2]))*constants[3]
    algebraic[33] = ((constants[0]*constants[1])/constants[2])*log(constants[10]/states[1])
    algebraic[44] = 0.100000/(1.00000+exp(0.0600000*((states[0]-algebraic[33])-200.000)))
    algebraic[45] = (3.00000*exp(0.000200000*((states[0]-algebraic[33])+100.000))+exp(0.100000*((states[0]-algebraic[33])-10.0000)))/(1.00000+exp(-0.500000*(states[0]-algebraic[33])))
    algebraic[46] = algebraic[44]/(algebraic[44]+algebraic[45])
    algebraic[47] = constants[13]*algebraic[46]*(power(constants[10]/5.40000, 1.0/2))*(states[0]-algebraic[33])
    algebraic[54] = constants[20]*states[16]*states[15]*(states[0]-algebraic[33])
    algebraic[48] = constants[14]*(power(constants[10]/5.40000, 1.0/2))*states[4]*states[5]*(states[0]-algebraic[33])
    algebraic[41] = ((constants[0]*constants[1])/constants[2])*log((constants[10]+constants[9]*constants[11])/(states[1]+constants[9]*states[2]))
    algebraic[49] = constants[15]*(power(states[6], 2.00000))*(states[0]-algebraic[41])
    algebraic[52] = (((constants[18]*states[11]*states[12]*states[13]*states[14]*4.00000*(states[0]-15.0000)*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(0.250000*states[10]*exp((2.00000*(states[0]-15.0000)*constants[2])/(constants[0]*constants[1]))-constants[12]))/(exp((2.00000*(states[0]-15.0000)*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[43] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[12]/states[3])
    algebraic[53] = constants[19]*(states[0]-algebraic[43])
    algebraic[58] = (constants[32]*(states[0]-algebraic[33]))/(1.00000+exp((25.0000-states[0])/5.98000))
    algebraic[57] = (constants[30]*states[3])/(states[3]+constants[31])
    algebraic[12] = stimulus
    rates[0] = (-1.00000/1.00000)*(algebraic[47]+algebraic[54]+algebraic[48]+algebraic[49]+algebraic[52]+algebraic[55]+algebraic[50]+algebraic[51]+algebraic[56]+algebraic[53]+algebraic[58]+algebraic[57]+algebraic[12])
    rates[1] = ((-1.00000*((algebraic[47]+algebraic[54]+algebraic[48]+algebraic[49]+algebraic[58]+algebraic[12])-2.00000*algebraic[55]))/(1.00000*constants[4]*constants[2]))*constants[3]
    algebraic[59] = constants[44]/(1.00000+(power(constants[42], 2.00000))/(power(states[3], 2.00000)))
    algebraic[60] = constants[43]*(states[17]-states[3])
    algebraic[61] = constants[41]*(states[10]-states[3])
    algebraic[63] = 1.00000/(1.00000+(constants[45]*constants[46])/(power(states[3]+constants[46], 2.00000)))
    rates[3] = algebraic[63]*((((algebraic[60]-algebraic[59])*constants[51])/constants[4]+algebraic[61])-(1.00000*((algebraic[53]+algebraic[57])-2.00000*algebraic[56])*constants[3])/(2.00000*1.00000*constants[4]*constants[2]))
    algebraic[62] = constants[38]-(constants[38]-constants[39])/(1.00000+power(constants[37]/states[17], 2.00000))
    algebraic[65] = constants[34]*algebraic[62]
    rates[18] = -algebraic[65]*states[10]*states[18]+constants[36]*(1.00000-states[18])
    algebraic[64] = constants[33]/algebraic[62]
    algebraic[66] = (algebraic[64]*(power(states[10], 2.00000))*states[18])/(constants[35]+algebraic[64]*(power(states[10], 2.00000)))
    algebraic[67] = constants[40]*algebraic[66]*(states[17]-states[10])
    algebraic[68] = 1.00000/(1.00000+(constants[47]*constants[48])/(power(states[17]+constants[48], 2.00000)))
    rates[17] = algebraic[68]*(algebraic[59]-(algebraic[67]+algebraic[60]))
    algebraic[69] = 1.00000/(1.00000+(constants[49]*constants[50])/(power(states[10]+constants[50], 2.00000)))
    rates[10] = algebraic[69]*(((-1.00000*algebraic[52]*constants[3])/(2.00000*1.00000*constants[52]*constants[2])+(algebraic[67]*constants[51])/constants[52])-(algebraic[61]*constants[4])/constants[52])
    return (rates, algebraic)
