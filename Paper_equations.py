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

def runge_kutta4(f, t, y, h, m, q_0):
    k1 = h * f(t, y, m, q_0)
    k2 = h * f(t + 0.5 * h, y + 0.5 * k1, m, q_0)
    k3 = h * f(t + 0.5 * h, y + 0.5 * k2, m, q_0)
    k4 = h * f(t + h, y + k3, m, q_0)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

def solve_with_odeint(f, t0, y0, tf, h, m, q_0):
    t_values = np.linspace(t0, tf, int(tf/h) + 1)  # Create time grid with more precision
    return odeint(f, y0, t_values, args=(m, q_0))

def paper_differential_equations(y, t, m, q_0):
    P_J = y[0]
    psi = y[1]
    J = y[2]
    theta = y[3]
    
    # hypergeometric function parameters definition
    q_wall = 3.5
    n = 2

    psi_p = psi/q_0*calc_hyp2f1(psi, q_0, q_wall, n)
    q = 1 #q_0 *(1 + (psi)**2)**(1/2)

    der_theta = - ((1-(2*psi)*np.cos(theta))*(P_J + psi_p))/np.sqrt(2*psi)*np.cos(theta) - m*np.cos(theta)/np.sqrt(2*psi) + (1-psi*np.cos(theta))**2*(P_J +psi_p)*q
    der_psi = -(P_J + psi_p)**2*(1 - np.sqrt(psi))*np.sqrt(psi)*np.sin(theta) - m*psi*np.sin(theta)
    der_P_J = 0
    der_J = (P_J + psi_p)*(1 - psi*np.cos(theta))

    return np.array([der_P_J, der_psi, der_J, der_theta])


def solve_runge_kutta4(f, t0, y0, tf, h, m, q_0):
    t_values = np.arange(t0, tf + h, h)
    y_values = [y0]
    t = t0
    y = y0
    for next_t in t_values[1:]:
        y = runge_kutta4(f, t, y, h, m, q_0)
        t = next_t
        y_values.append(y)
    return np.array(t_values), np.array(y_values)
