# This file defines the controllers ODE solutions to be imported elsewhere

def dol_solution(params):
    """
    function to get data frame for open loop
    """
    alpha_X      = params.get("alpha_X", 0.1)
    alpha_Y      = params.get("alpha_Y", 0.001)
    alpha_Y_plus = params.get("alpha_Y+", 1.0)    
    delta_X      = params.get("delta_X", 0.0005)
    delta_Y      = params.get("delta_Y", 0.0005)
    kappa        = params.get("kappa", 0.0)
    omega        = params.get("omega", 5e6)
    v            = params.get("nu", 1.0)
    alpha_O      = params.get("alpha_O", 0.001)
    alpha_O_plus = params.get("alpha_O+", 1.0)    
    delta_O      = params.get("delta_O", 0.002)

    # Initial conditions
    P_y = 20e-9
    X0  = 0
    Y0  = 0
    Pp0 = 0
    O0  = 0
    variable0  = [X0, Y0, Pp0, O0]
    
    # Time vector
    t_min   = np.linspace(0, 600, 1000)
    t_sec   = t_min * 60.0
    
    def direct_open_loop(variable, t):
        X, Y, P_plus, O = variable
        if t >= 3600 :
            P_x = 10e-9
        else:
            P_x = 0
            
        dX_dt = alpha_X * P_x - delta_X * X - kappa * X * Y - omega * X * P_y + v * P_plus
        dY_dt = alpha_Y * P_y + alpha_Y_plus * P_plus - delta_Y * Y - kappa * X * Y
        dPp_dt = omega * X * P_y - v * P_plus
        dO_dt = alpha_O * P_y + alpha_O_plus * P_plus - delta_O * O
        
        return [dX_dt, dY_dt, dPp_dt, dO_dt]
    
    # Solve ODE system
    
    sol = odeint(direct_open_loop, variable0, t_sec)
    X, Y, P_plus, O = sol.T
    
    # Build and return DataFrame
    df = pd.DataFrame({
        "time[s]": t_sec,
        "X":       X,
        "Y":       Y,
        "P_Y+":    P_plus,
        "O":       O
    })
    return df

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

def iol_solution():
    return

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


