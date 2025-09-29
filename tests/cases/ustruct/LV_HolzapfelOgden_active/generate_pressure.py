import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

params = {}

params['alpha_min'] = -30
params['alpha_max'] = 5
params['alpha_pre'] = 5
params['alpha_mid'] = 1
params['sigma_pre'] = 7000
params['sigma_mid'] = 16000
params['sigma_0'] = 1.5e5
params['t_sys_pre'] = 0.17
params['t_dias_pre'] = 0.484
params['t_sys'] = 0.16
params['t_dias'] = 0.484
params['gamma'] = 0.005

def S_plus(delta_t, params):
    return 0.5* (1 + np.tanh(delta_t/params['gamma']))

def S_minus(delta_t, params):
    return 0.5* (1 - np.tanh(delta_t/params['gamma']))

def g_pre(t, params):
    return S_minus(t - params['t_dias_pre'], params)

def f(t, params):
    return S_plus(t - params['t_sys'], params) * S_minus(t - params['t_dias'], params)

def f_pre(t, params):
    return S_plus(t - params['t_sys_pre'], params) * S_minus(t - params['t_dias_pre'], params)

def a(t, params):
    return params['alpha_max'] * f(t, params) + params['alpha_min'] * (1-f(t, params))

def a_pre(t, params):
    return params['alpha_max'] * f_pre(t, params) + params['alpha_min'] * (1-f_pre(t, params))

def b(t, params):
    return a_pre(t, params) + params['alpha_pre'] * g_pre(t, params) + params['alpha_mid']

# Differential equation for pressure
def dpdt(t, p):
    b_abs = abs(b(t, params))
    b_max = max(b(t, params), 0)
    g_pre_max = max(g_pre(t, params), 0)
    return -b_abs*p + params['sigma_mid']*b_max + params['sigma_pre']*g_pre_max

# Differential equation for stress
def dtaudt(t, tau):
    a_abs = abs(a(t, params))
    a_max = max(a(t, params), 0)
    return -a_abs*tau + params['sigma_0']*a_max


def save_sol(sol, filename, f_modes = 512):
    # save data to .dat file
    # first row is [len(sol.y[0]) f_modes]
    # first column is sol.t rounded to six decimal places
    # second column is sol.y[0] rounded to six decimal places

    data = np.column_stack((np.round(sol.t, 6), np.round(sol.y[0], 6)))

    data = np.vstack((np.array([len(sol.y[0]), f_modes]), data))

    # Save to file with the first row as integers and the rest as floats
    with open(filename, 'w') as f:
        # Write the first row as integers
        np.savetxt(f, [data[0]], fmt='%d', delimiter=' ')
        # Write the remaining rows as floats with 6 decimal places
        np.savetxt(f, data[1:], fmt='%1.6f', delimiter=' ')

# Solve the ODE for pressure

p0 = 0
t_eval = np.linspace(0, 1, 1001)
p_sol = solve_ivp(dpdt, [0, 1], [p0], t_eval=t_eval, method='DOP853')
print('Max pressure:', max(p_sol.y[0]))

plt.plot(p_sol.t, p_sol.y[0])
plt.xlabel('Time')
plt.ylabel('Pressure')
'''
# load pressure_original.dat and convert string to float
p_original = np.loadtxt('pressure_original.dat', delimiter=',')
p_original[1:,1] = 0.1*p_original[1:,1]
plt.plot(p_original[1:,0], p_original[1:,1])

with open('pressure_original_Pa.dat', 'w') as f:
        # Write the first row as integers
        np.savetxt(f, [p_original[0]], fmt='%d', delimiter=' ')
        # Write the remaining rows as floats with 6 decimal places
        np.savetxt(f, p_original[1:], fmt='%1.6f', delimiter=' ')
'''
plt.show()


# Save the solution to a file
save_sol(p_sol, 'pressure.dat', 512)

# Solve the ODE for stress

tau0 = 0
t_eval = np.linspace(0, 1, 1001)
tau_sol = solve_ivp(dtaudt, [0, 1], [tau0], t_eval=t_eval)

print('Max stress:', max(tau_sol.y[0]))

plt.plot(tau_sol.t, tau_sol.y[0])
plt.xlabel('Time')
plt.ylabel('Stress')
'''
# load stress_original.dat and convert string to float
tau_original = np.loadtxt('stress_original.dat', delimiter=' ')
tau_original[1:,1] = 0.1*tau_original[1:,1]
plt.plot(tau_original[1:,0], tau_original[1:,1])

with open('stress_original_Pa.dat', 'w') as f:
        # Write the first row as integers
        np.savetxt(f, [tau_original[0]], fmt='%d', delimiter=' ')
        # Write the remaining rows as floats with 6 decimal places
        np.savetxt(f, tau_original[1:], fmt='%1.6f', delimiter=' ')
'''
plt.show()

# Save the solution to a file
save_sol(tau_sol, 'LV_stress.dat', 512)


