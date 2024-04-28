# Generalized Coordinates and Hamiltonian Dynamics

This repository contains Python code to solve ordinary differential equations (ODEs) describing the equations of motion in the GC of a particle. The equations are solved numerically using the `odeint` function from the SciPy library.

## Equations

The equations are derived from the Hamiltonian of the GC. The equations used in this code are:

```math
\begin{align*}
\frac{dP_{\zeta}}{dt} &= 0 \\
\frac{d\psi}{dt} &= -\left(P_{\zeta} + \frac{\psi}{q_0}F\left(\frac{\psi}{q_0}\right)\right)^2\left(1 - \sqrt{2\psi}\cos(\theta)\right)\sqrt{2\psi}\sin(\theta) - m\sqrt{2\psi}\sin(\theta) \\
\frac{d\zeta}{dt} &= (P_{\zeta} + \frac{\psi}{q_0}F\left(\frac{\psi}{q_0}\right))\left(1 - \sqrt{2\psi}\cos(\theta)\right)^2 \\
\frac{d\theta}{dt} &= -\frac{\cos(\theta)}{\sqrt{2\psi}}(1 - \sqrt{2\psi}\cos(\theta))\left(P_{\zeta} + \frac{\psi}{q_0}F\left(\frac{\psi}{q_0}\right)\right)^2 - \frac{m\cos(\theta)}{{\sqrt{2\psi}}} + \frac{(1 - \sqrt{2\psi}\cos(\theta))^2(P_{\zeta} + \frac{\psi}{q_0}F\left(\frac{\psi}{q_0}\right))}{q}
\end{align*}
```

Where:
- $P_{\zeta}$ is the toroidal momentum.
- $\psi = P_{\theta}$ in the LAR approximation represents the poloidal momentum.
- $\zeta$ is the toroidal angle.
- $\theta$ is poloidal angle.
- $\mu$ is the magnetic moment.
- $q_0$ is a constant.
- $F$ is the hypergeometric function.

## Hamiltonian Equation

The Hamiltonian $W$ with respect to $(\theta, \zeta, \psi, P_{\zeta})$ coordinates for the system is given by:

```math
W = \left(P_J + \frac{\psi}{q_0}F\left(\frac{\psi}{q_0}\right)\right)^2 \left(1 - \sqrt{2\psi}\cos(\theta)\right)^2 / 2 + m\left(1 - \sqrt{2\psi}\cos(\theta)\right)
```
