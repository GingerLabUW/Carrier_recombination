import numpy as np
from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from scipy.special import gamma

def generate_single_exp(t, tau):
    I = np.exp(-t/tau)
    return I

def fun(y, t, k1, k2=8e-11):
    #y1, y2 = y #y1 is the carrier concentration and y2 is the photon concentration in the film
    dy1dt = -k1*y - k2*y**2 #+ (c/n)*Alpha_avg*y2 #change in the carrier population as function of time
    #dy2dt = -(c/n)*Alpha_avg*y2 + (k_b*y1**2*P_stay) #change in the photon density as function of time
    return dy1dt

# Solve the coupled ODE's
def ode_solve(N0, k1, k2=8e-11):#, k_A):
    sim, infod = odeint(fun, N0 , t, args=(k1, k2), full_output = 1, mxstep=5000000)
    return sim

def stretch_exp_fit(TRPL, t, Tc = (0,1e4*1e-9), Beta = (0,1), A = (0,1.5)):

    def exp_stretch(t, tc, beta, a):
        return a * np.exp(-((1.0 / tc) * t) ** beta)

    def avg_tau_from_exp_stretch(tc, beta):
        return (tc / beta) * gamma(1.0 / beta)

    def Diff_Ev_Fit_SE(TRPL):

        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tc = params[0]
            beta = params[1]
            a = params[2]
            
            PL_sim = exp_stretch(t,tc,beta,a)

            Resid= (np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data

        bounds = [Tc, Beta, A]

        result = differential_evolution(residuals, bounds)
        return result.x

    p = Diff_Ev_Fit_SE(TRPL)

    tc = p[0]
    beta = p[1]
    a = p[2]
    
    PL_fit = exp_stretch(t,tc,beta,a)

    avg_tau = avg_tau_from_exp_stretch(tc,beta)

    return tc, beta, a, avg_tau, PL_fit

def double_exp_fit(TRPL, t, tau1_bounds=(0,1000*1e-9), a1_bounds=(0,1), tau2_bounds=(0,10000*1e-9), a2_bounds=(0,1)):

    def single_exp(t, tau, a):
        return (a * np.exp(-((1.0 / tau)*t)))

    def double_exp(t, tau1, a1, tau2, a2):
        return ((single_exp(t, tau1, a1)) + (single_exp(t, tau2, a2)))

    def avg_tau_from_double_exp(tau1, a1, tau2, a2):
        return (((tau1*a1) + (tau2*a2))/(a1+a2))

    def Diff_Ev_Fit_DE(TRPL):

        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tau1 = params[0]
            a1 = params[1]
            tau2 = params[2]
            a2 = params[3]


            PL_sim = double_exp(t,tau1, a1, tau2, a2)

            Resid= (np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data

        bounds = [tau1_bounds, a1_bounds, tau2_bounds, a2_bounds]

        result = differential_evolution(residuals, bounds)
        return result.x

    p = Diff_Ev_Fit_DE(TRPL)

    tau1 = p[0]
    a1 = p[1]
    tau2 = p[2]
    a2 = p[3]

    
    PL_fit = double_exp(t, tau1, a1, tau2, a2)

    avg_tau = avg_tau_from_double_exp(tau1, a1, tau2, a2)

    return tau1, a1, tau2, a2, avg_tau, PL_fit

def single_exp_fit(TRPL, t, tau_bounds=(0,10000*1e-9), a_bounds=(0,1)):

    def single_exp(t, tau, a):
        return (a * np.exp(-((1.0 / tau)*t)))

    def avg_tau_from_single_exp(tau, a):
        return ((tau*a)/(a))

    def Diff_Ev_Fit_SiE(TRPL):

        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tau = params[0]
            a = params[1]

            PL_sim = single_exp(t,tau, a)

            Resid= (np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data

        bounds = [tau_bounds, a_bounds]

        result = differential_evolution(residuals, bounds)
        return result.x

    p = Diff_Ev_Fit_SiE(TRPL)

    tau = p[0]
    a = p[1]



    PL_fit = single_exp(t, tau, a)

    avg_tau = avg_tau_from_single_exp(tau, a)

    return tau, a, avg_tau, PL_fit

def calculate_surface_lifetime(effective_lifetime, bulk_lifetime=8000):
    """Calculate surface lifetime for a semiconductor given effective lifetime and bulk lifetime

    effective_lifetime : preferably in nanoseconds, units can be different from nanoseconds but bulk_lifetime must be modified accordingly
    bulk_lifetime : 8000, in nanoseconds, can be modified
    """
    surface_lifetime = (effective_lifetime * bulk_lifetime)/(bulk_lifetime - effective_lifetime)
    return surface_lifetime

def calculate_srv (surface_lifetime, diffusion_coefficient = 0.9, thickness = 400):
    """
    Calculate SRV for two different conditions (SRV1=0 and SRV1=SRV2) given surface lifetime (in ns), diffusion coefficient (in cm2/s) and thickness (in nm)

    Returns : SRV for SRV1=0 condition, SRV for SRV1=SRV2 condition
    """

    thickness = thickness*1e-7 # convert to cm
    diffusion_coefficient = diffusion_coefficient # in cm2/s

    srv1_srv2_equal = thickness / (2*((1e-9*surface_lifetime) - ((1/diffusion_coefficient)*((thickness/np.pi)**2)) ))
    srv1_zero = thickness / ((1e-9*surface_lifetime) - ((4/diffusion_coefficient)*((thickness/np.pi)**2)) )

    return srv1_zero, srv1_srv2_equal

def triple_exp_fit(TRPL, t, tau1_bounds=(0,1000*1e-9), a1_bounds=(0,1), tau2_bounds=(0,10000*1e-9), a2_bounds=(0,1), tau3_bounds=(0,10000*1e-9), a3_bounds=(0,1)):

    def single_exp(t, tau, a):
        return (a * np.exp(-((1.0 / tau)*t) ))

    def avg_tau_from_triple_exp(tau1, a1, tau2, a2, tau3, a3):
        return ((tau1*a1) + (tau2*a2) + (tau3*a3))/(a1+a2+a3)

    def triple_exp(t, tau1, a1, tau2, a2, tau3, a3):
        return ((single_exp(t, tau1, a1)) + (single_exp(t, tau2, a2)) + (single_exp(t, tau3, a3)))

    def Diff_Ev_Fit_TE(TRPL):

        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tau1 = params[0]
            a1 = params[1]
            tau2 = params[2]
            a2 = params[3]
            tau3= params[4]
            a3 = params[5]
         

            PL_sim = triple_exp(t,tau1, a1, tau2, a2, tau3, a3)

            Resid= (np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data

        bounds = [tau1_bounds, a1_bounds, tau2_bounds, a2_bounds, tau3_bounds, a3_bounds]

        result = differential_evolution(residuals, bounds)
        return result.x

    p = Diff_Ev_Fit_TE(TRPL)

    tau1 = p[0]
    a1 = p[1]
    tau2 = p[2]
    a2 = p[3]
    tau3= p[4]
    a3 = p[5]
    


    PL_fit = triple_exp(t, tau1, a1, tau2, a2, tau3, a3)

    avg_tau = avg_tau_from_triple_exp(tau1, a1, tau2, a2, tau3, a3)

    return tau1, a1, tau2, a2, tau3, a3, avg_tau, PL_fit

def dynamic_trap_fit(TRPL, t, Tc = (0,1e4*1e-9), Beta = (0,1), A = (0,100), tau_bounds=(0,1e4*1e-9), a_bounds=(0,100)):

    def dynamic_trap(t, tc, beta, a1, tau, a2):
        return (a1 * np.exp(-((1.0 / tc) * t) ** beta)) + (a2 * np.exp(-((1.0 / tau)*t) ))

    def avg_tau_from_exp_stretch(tc, beta):
        return (tc / beta) * gamma(1.0 / beta)

    def Diff_Ev_Fit_DT(TRPL):

        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tc = params[0]
            beta = params[1]
            a1 = params[2]
            tau = params[3]
            a2 = params[4]
                
            PL_sim = dynamic_trap(t,tc,beta,a1,tau,a2)

            Resid= (np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data

        bounds = [Tc, Beta, A, tau_bounds, a_bounds]

        result = differential_evolution(residuals, bounds)
        return result.x

    p = Diff_Ev_Fit_DT(TRPL)

    tc = p[0]
    beta = p[1]
    a1 = p[2]
    tau = p[3]
    a2= p[4]
                
    
    PL_fit = dynamic_trap(t,tc,beta,a1,tau,a2)

    avg_tau = avg_tau_from_exp_stretch(tc,beta)

    return tc, beta, a1, tau, a2, avg_tau, PL_fit
