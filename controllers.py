# This file defines the controllers ODE solutions to be imported elsewhere

def dol_solution(conds):
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="direct_closed_loop", temperature_C=temp)
    
    t = np.linspace(0, conds["t_end"], conds["t_steps"])
    
    def open_system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"] - parameters["kappa"]*X*Y - parameters["omega"]*X*PY + parameters["nu"]*PY_plus
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = 0
        dPY_plusdt = parameters["omega"]*X*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    
    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
    })
    
    return response

def dcl_solution(conds):
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="direct_closed_loop", temperature_C=temp)
    
    t = np.linspace(0, conds["t_end"], conds["t_steps"])
    
    def system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"] - parameters["kappa"]*X*Y - parameters["omega"]*X*PY + parameters["nu"]*PY_plus
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = 0
        dPY_plusdt = parameters["omega"]*X*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    
    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
    })
    
    return response

def iol_solution(conds):
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="direct_closed_loop", temperature_C=temp)
    
    t = np.linspace(0, conds["t_end"], conds["t_steps"])
    
    def system(vars, t):
        X, Y, A, PY_plus, O = vars

        dXdt = parameters["alpha_X"]*PX - parameters["delta_X"] *X - parameters["kappa"]*X*Y
        dYdt = parameters["alpha_Y"]*PY + parameters["alpha_Y+"]*PY_plus - parameters["delta_Y"]*Y - parameters["kappa"]*X*Y
        dAdt = parameters["beta_A"] *X - parameters["delta_A"] *A - parameters["omega"] *A *PY + parameters["nu"] * PY_plus
        dPY_plusdt = parameters["omega"]*A*PY - parameters["nu"]*PY_plus
        dOdt = parameters["alpha_O"]*PY + parameters["alpha_O+"]*PY_plus - parameters["delta_O"]*O
        return [dXdt, dYdt, dAdt, dPY_plusdt, dOdt]
    
    initial_conditions = [X0, Y0, A0, PY_plus0, O0]
    
    solution = scipy.integrate.odeint(system, initial_conditions, t)
    
    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
    })
    
    return response

def icl_solution(conds):
    temp = conds["temperature"]
    PX = conds["PX"]
    PY = conds["PY"]
    
    X0 = conds["X_init"]
    Y0 = conds["Y_init"]
    A0 = conds["A_init"]
    PY_plus0 = conds["PY_plus_init"]
    O0 = conds["O_init"]
    
    parameters = rate(controller="indirect_closed_loop", temperature_C=temp)
    
    t = np.linspace(0, conds["t_end"], conds["t_steps"])
    
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
    
    response = pd.DataFrame({
        "time[s]": t,
        "X": solution[:, 0],
        "Y": solution[:, 1],
        "A": solution[:, 2],
        "PY_+": solution[:, 3],
        "O": solution[:, 4],
    })
    return response


