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

def icl_solution():
    return


def masterEQ(xy, t,temp):
    X, Y, A, P_Y_plus, O, P_X,P_Y = xy
    rate_dict = rate(controller="indirect_closed_loop", temperature_C=temp)
    dX_dt = (rate_dict["alpha_X"] * P_X) - (rate_dict["delta_X"] * X) - (rate_dict["kappa"] * X * Y)
    dY_dt = (rate_dict["alpha_Y"] * P_Y) + (rate_dict["alpha_Y+"] * P_Y_plus) - (rate_dict["delta_Y"] * Y) - (rate_dict["kappa"] * X * Y)
    dA_dt = (rate_dict["beta_A"] * X) - (rate_dict["delta_A"] * A) - (rate_dict["omega"] * A * P_Y) + (rate_dict["nu"] * P_Y_plus)
    dP_Y_plus = (rate_dict["omega"] * A * P_Y) - (rate_dict["nu"] * P_Y_plus)
    dO_dt = (rate_dict["alpha_O"] * P_Y) + (rate_dict["alpha_O+"] * P_Y_plus) - (rate_dict["delta_O"] * O)
    dP_X_dt = 0
    dP_Y_dt = 0
    return np.array([dX_dt, dY_dt, dA_dt, dP_Y_plus, dO_dt, dP_X_dt, dP_Y_dt])

def indirect_closed(eq, params):
    t = np.linspace(0, 100, 1000)
    xyapoXY_0 = np.array([0,0,0,0,0,10,20])
    teto = odeint(eq, xyapoXY_0, t,  args = (params,))
    return teto