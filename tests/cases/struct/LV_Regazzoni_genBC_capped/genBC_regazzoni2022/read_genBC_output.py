
def read_genBC_output(output_file, parameters_file):
    import numpy as np

    # Read 0D-0D data in data
    data = np.loadtxt(output_file)

    # Read data into dictionary
    results_dict = {
        'p_LV': data[:, 0],  # mmHg
        'p_LV_cap': data[:, 1],
        'V_RV': data[:, 2],
        'V_LA': data[:, 3],  # mL
        'V_RA': data[:, 4],
        'p_AR_SYS': data[:, 5],
        'p_VEN_SYS': data[:, 6],
        'p_AR_PUL': data[:, 7],
        'p_VEN_PUL': data[:, 8],
        'Q_AR_SYS': data[:, 9],  # mL/s
        'Q_VEN_SYS': data[:, 10],
        'Q_AR_PUL': data[:, 11],
        'Q_VEN_PUL': data[:, 12],
        'V_LV_0D': data[:, 13],
        'time': data[:, 14],
        'R_MV': data[:, 15],
        'R_AV': data[:, 16],
        'R_TV': data[:, 17],
        'R_PV': data[:, 18],
        'Q_LV': data[:, 19],
        'A_RV': data[:, 20],
        'A_LA': data[:, 21],
        'A_RA': data[:, 22]
    }

    # Read parameters
    parameters = {}
    with open(parameters_file, 'r') as file:
        for line in file:
            # Remove comments
            line = line.split('!')[0]
            if line.startswith('      REAL(KIND=8), PARAMETER ::'):
                parts = line.split('::')
                var_name = parts[1].split('=')[0].strip()
                try:
                    var_value = eval(parts[1].split('=')[1].strip())
                    parameters[var_name] = var_value
                except:
                    # If eval fails, check if the parts of the right-hand side are keys in the parameters dictionary
                    rhs_parts = parts[1].split('=')[1].strip().split()
                    for i, part in enumerate(rhs_parts):
                        if part in parameters:
                            rhs_parts[i] = str(parameters[part])
                    # Try to evaluate the right-hand side again
                    try:
                        var_value = eval(' '.join(rhs_parts))
                        parameters[var_name] = var_value
                    except:
                        print('Could not read parameters from line: {}'.format(line))

    # Add parameters to results_dict
    results_dict['parameters'] = parameters

    # Compute additional quantities
    # RV and Atrial pressures
    V0_RV = parameters['V0_RV']  # mL
    V0_LA = parameters['V0_LA']  # mL
    V0_RA = parameters['V0_RA']  # mL
    E_RV_pas = parameters['E_RV_pas']  # mmHg/mL
    E_RV_act = parameters['E_RV_act']  # mmHg/mL
    E_LA_pas = parameters['E_LA_pas']  # mmHg/mL
    E_LA_act = parameters['E_LA_act']  # mmHg/mL
    E_RA_pas = parameters['E_RA_pas']  # mmHg/mL
    E_RA_act = parameters['E_RA_act']  # mmHg/mL

    results_dict['p_RV'] = (E_RV_pas + E_RV_act * results_dict['A_RV']) * (results_dict['V_RV'] - V0_RV)  # mmHg
    results_dict['p_LA'] = (E_LA_pas + E_LA_act * results_dict['A_LA']) * (results_dict['V_LA'] - V0_LA)  # mmHg
    results_dict['p_RA'] = (E_RA_pas + E_RA_act * results_dict['A_RA']) * (results_dict['V_RA'] - V0_RA)  # mmHg

    # Valve flow rates
    results_dict['Q_MV'] = (results_dict['p_LA'] - results_dict['p_LV']) / results_dict['R_MV']  # Mitral valve forward flow
    results_dict['Q_AV'] = (results_dict['p_LV'] - results_dict['p_AR_SYS']) / results_dict['R_AV']  # Aortic valve forward flow
    results_dict['Q_TV'] = (results_dict['p_RA'] - results_dict['p_RV']) / results_dict['R_TV']  # Tricuspid valve forward flow
    results_dict['Q_PV'] = (results_dict['p_RV'] - results_dict['p_AR_PUL']) / results_dict['R_PV']  # Pulmonary valve forward flow

    # Compute vascular volumes from pressure solution
    results_dict['V_AR_SYS'] = results_dict['p_AR_SYS'] * parameters['C_AR_SYS']
    results_dict['V_VEN_SYS'] = results_dict['p_VEN_SYS'] * parameters['C_VEN_SYS']
    results_dict['V_AR_PUL'] = results_dict['p_AR_PUL'] * parameters['C_AR_PUL']
    results_dict['V_VEN_PUL'] = results_dict['p_VEN_PUL'] * parameters['C_VEN_PUL']

    # Compute total blood volume
    results_dict['V_tot'] = results_dict['V_LA'] + results_dict['V_RA'] + results_dict['V_LV_0D'] + results_dict['V_RV'] + \
                           results_dict['V_AR_SYS'] + results_dict['V_VEN_SYS'] + results_dict['V_AR_PUL'] + results_dict['V_VEN_PUL']


    return results_dict

if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(__file__))
    output_file = '../AllData'
    parameters_file = 'include/parameters.f'
    results_dict = read_genBC_output(output_file, parameters_file)
    print(results_dict)