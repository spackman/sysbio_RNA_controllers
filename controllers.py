import numpy as np
import pandas as pd
import scipy
from math import exp
from typing import Dict, List, Tuple


# ==========================================================
# UTILITY FUNCTIONS
# ==========================================================

def arrhenius_multiplication_factor(temperature_C : List , activation_energy_kcal : float ) -> List:
    temp_K = temperature_C + 273.15 # [ K ]
    reference_temp = 298.15 # [ K ]
    gas_constant = 1.9872E-3 # [ kcal / mol-K]
    factor = exp(activation_energy_kcal*(temp_K - reference_temp) / (gas_constant * reference_temp * temp_K)) # [ unitless ]
    return factor

def rate(controller: str, temperature_C: float = 25) -> Dict[str, float]:
    """
    Get the dictionary of model parameters for the selected controller type
    at a given temperature (in Celsius). Parameters are scaled using the
    Arrhenius equation based on their reaction type (enzyme or RNA).

    Parameters:
    - controller (str): The controller type, must be one of:
        ["direct_open_loop", "direct_closed_loop", "indirect_open_loop", "indirect_closed_loop"]
    - temperature_C (float): Temperature in degrees Celsius at which to scale the parameters.
                             Default is 25C, which references the parameters used in the paper.

    Returns:
    - Dict[str, float]: 
        - A dictionary of parameter names and their temperature-scaled values.
        
    Raises:
    - TypeError: If the controller type is not recognized.
    """
    controllers = ["direct_open_loop",
                   "direct_closed_loop",
                   "indirect_open_loop",
                   "indirect_closed_loop"]
    
    if controller not in controllers:
        print("The selected controller is not supported. Select a controller from the list:\n")
        print("\n".join(controllers))
        raise TypeError(f"Controller '{controller}' not supported.")

    parameters = pd.read_csv("kinetic_parameters.csv")

    enzyme_rxn_Ea = 15 # [ kcal / mol]
    RNA_rxn_Ea = 25 # [ kcal / mol ]

    enzyme_factor = arrhenius_multiplication_factor(temperature_C, activation_energy_kcal=enzyme_rxn_Ea)
    RNA_factor = arrhenius_multiplication_factor(temperature_C, activation_energy_kcal=RNA_rxn_Ea)

    def scale_row(row):
        factor = enzyme_factor if row["type"] == "enzyme" else RNA_factor
        value = row[controller] * factor
        return float(f"{value:.5g}") # round to 5 sig figs

    parameters["scaled_value"] = parameters.apply(scale_row, axis=1)

    rate_at_temperature = dict(zip(parameters["parameter"], parameters["scaled_value"]))

    for param in ["kappa", "omega"]:
        # convert concentration dependent parameters to  1/ nM-s from 1 / M-s
        # so that all concentrations can be in nM
        rate_at_temperature[param] = rate_at_temperature[param] / 1E9
    return rate_at_temperature
    

def to_latex_label(param_name: str) -> str:
    """
    Convert a parameter name like 'alpha_X' or 'kappa' into a LaTeX-formatted
    label for matplotlib plotting. 
    Edited from GPT prompt: 
    'write a small function that converts strings to properly formatted 
    latex to be displayed in matplotlib like r"\alpha_X" and "\omega" etc'

    Examples:
    - 'alpha_X' -> r"$\\alpha_X$"
    - 'kappa'   -> r"$\\kappa$"
    """
    greek_letters = {'alpha', 'beta', 'delta','kappa', 'nu', 'omega' }

    parts = param_name.split('_')

    if parts[0] in greek_letters:
        base = f"\\{parts[0]}"
    else:
        base = parts[0]  # Not a greek letter, leave as is

    if len(parts) == 1:
        return rf"${base}$"
    elif len(parts) == 2:
        return rf"${base}_{{{parts[1]}}}$"
    else:
        # Join all remaining parts with _ and allow for sub/subscripts
        subscript = "_".join(parts[1:])
        return rf"${base}_{{{subscript}}}$"
    

# ==========================================================
# ODE SOLUTIONS
# ==========================================================

def dol_solution(conds):
    ''' Direct Open Loop Controller'''
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="direct_open_loop", temperature_C=temp)
    
    t = np.linspace(0, conds["t_end"], conds["t_steps"])
    
    def system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"]*X - parameters["kappa"]*X*Y - parameters["omega"]*X*PY + parameters["nu"]*PY_plus
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = 0
        dPY_plusdt = parameters["omega"]*X*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    reference = parameters["alpha_O+"] * parameters["alpha_X"] / (parameters["delta_O"] * parameters["alpha_Y+"]) * PX
    
    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
        "reference": reference
    })
    
    return response

def dcl_solution(conds):
    ''' Direct Closed Loop Controller'''
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="direct_closed_loop", temperature_C=temp)
    
    t = np.linspace(0, int(conds["t_end"]), int(conds["t_steps"]))
    
    def system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"]*X - parameters["kappa"]*X*Y - parameters["omega"]*X*PY + parameters["nu"]*PY_plus
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = 0
        dPY_plusdt = parameters["omega"]*X*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    reference = parameters["alpha_O+"] * parameters["alpha_X"] / (parameters["delta_O"] * parameters["alpha_Y+"]) * PX

    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
        "reference": reference

    })
    
    return response

def iol_solution(conds):
    ''' Indirect Open Loop Controller'''
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="indirect_open_loop", temperature_C=temp)
    
    t = np.linspace(0, conds["t_end"], conds["t_steps"])
    
    def system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"]*X - parameters["kappa"]*X*Y
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = parameters["beta_A"]*X - parameters["delta_A"]*A - parameters["omega"]*A*PY + parameters["nu"]*PY_plus
        dPY_plusdt = parameters["omega"]*A*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    reference = parameters["alpha_O+"] * parameters["alpha_X"] / (parameters["delta_O"] * parameters["alpha_Y+"]) * PX

    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
        "reference" : reference
    })
    
    return response

def icl_solution(conds):
    ''' Indirect Closed Loop Controller'''
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="indirect_closed_loop", temperature_C=temp)
    
    t = np.linspace(0, int(conds["t_end"]), int(conds["t_steps"]))
    
    def system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"]*X - parameters["kappa"]*X*Y 
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = (parameters["beta_A"] * X) - (parameters["delta_A"] * A) - (parameters["omega"] * A * PY) + (parameters["nu"] *PY_plus)
        dPY_plusdt = parameters["omega"]*A*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    reference = parameters["alpha_O+"] * parameters["alpha_X"] / (parameters["delta_O"] * parameters["alpha_Y+"]) * PX

    
    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
        "reference" : reference
    })
    return response


