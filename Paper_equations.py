from scipy.special import hyp2f1
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


def calc_hyp2f1(psi, q_0, q_wall, n): 
    a = 1/n
    b = 1/n
    c = 1 + 1/n
    d = (1 - (q_wall/q_0)**n) * (psi)**n
    return hyp2f1(a, b, c, d)

def solve_with_odeint(f, t0, y0, tf, h, m, q_wall, q_0, n):
    t_values = np.linspace(t0, tf, int(tf/h) + 1)  
    return odeint(f, y0, t_values, args=(m, q_wall, q_0, n))

def paper_differential_equations(y, t, m, q_wall, q_0, n):
    # y      : solution for the ode 
    # t      : time parameter
    # m      : magnetic moment of order 10^-5
    # q_wall : 
    # q_0    : we consider it 1

    P_J = y[0]
    psi = y[1]
    J = y[2]
    theta = y[3]

    psi_p = psi/q_0*calc_hyp2f1(psi, q_0, q_wall, n)
    q = q_0*(1 + (psi)**2)**(1/2)

    der_theta = - np.cos(theta)*(1-np.sqrt(2*psi)*np.cos(theta))*((P_J + psi_p))**2/np.sqrt(2*psi) - m*np.cos(theta)/np.sqrt(2*psi) + (1 - np.sqrt(2*psi)*np.cos(theta))**2*(P_J +psi_p)/q
    der_psi = -(P_J + psi_p)**2*(1 - np.sqrt(2*psi)*np.cos(theta))*np.sqrt(2*psi)*np.sin(theta) - m*np.sqrt(2*psi)*np.sin(theta)
    der_P_J = 0
    der_J = (P_J + psi_p)*(1 - np.sqrt(2*psi)*np.cos(theta))**2

    return np.array([der_P_J, der_psi, der_J, der_theta])

