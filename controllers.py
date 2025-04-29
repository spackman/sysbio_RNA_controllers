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

def dcl_solution():
    return

def iol_solution():
    return

def icl_solution():
    return