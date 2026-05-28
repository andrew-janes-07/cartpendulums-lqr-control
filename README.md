# Cart-Pole LQR Control: Single, Double, and Triple Inverted Pendulums

Stabilizing pendulums in various unstable equilibriums on a cart via LQR with deadband, validated against
energy conservation to [Verify]**< 1e-9** in the uncontrolled system.

<p align="center">
  <img src="media/triple_stabilization.gif" width="600" alt="Triple pendulum stabilization">
</p>

## Overview

This project implements full nonlinear dynamics, linearization, and LQR-with-deadband
control for cart-pole systems with one, two, and three links. The triple-link case
is the headline result: an open-loop unstable 8-dimensional system with strong
coupling between links, stabilized at the upright equilibrium (and others) from non-trivial
initial perturbations.

**What's interesting here:**

- The Coriolis matrix for the *n*-link case is non-trivial to get right —
  energy conservation in the uncontrolled simulation is the test that catches
  almost every implementation bug.
- LQR alone produces high-frequency chatter near the equilibrium due to small
  state errors driving non-zero control. A deadband on the state error
  eliminates this chatter at the cost of a tiny steady-state offset, which is
  a reasonable trade for any real actuator.
- The triple pendulum is genuinely hard: the unstable manifold is
  three-dimensional and the deviations allowed before the system diverges from stability are small. 

## Results

### Stabilization from perturbation

[VERIFY: describe the perturbation magnitude you tested, e.g. "All three links
perturbed 15° from upright; cart starts at rest"]

<p align="center">
  <img src="media/triple_state_trajectories.png" width="700" alt="State trajectories">
</p>

All four positions (cart, θ₁, θ₂, θ₃) converge to the upright reference
(x = 0, θᵢ = π) within [VERIFY: ~3 seconds].

### Energy conservation and state space curves (uncontrolled, validation)

<p align="center">
  <img src="media/energy_conservation.png" width="600" alt="Energy drift">
</p>

Total mechanical energy drifts by less than [Verify 1e-9] over [VERIFY: 10 seconds] of
chaotic uncontrolled motion using `solve_ivp` with DOP853, rtol=1e-10,
atol=1e-12. This is the primary check that the inertia and Coriolis terms
are correct.

<p align="center">
  <img src="media/small_angle_coord_velocity_phase_space.png
</p>

The phase space plots of the small angle tests of the uncontrolled triple pendulum on a cart run for a t_span (0,1000)
trace out the expected regions the state space is allowed to reach. 


### Control effort

<p align="center">
  <img src="media/control_effort.png" width="600" alt="Control input over time">
</p>

The control function is continuous but only activates once above the deadband_on parameter. 
It then linear decreases to zero to the deadband_off parameter value if it dips below deadband_on. 
This approximates the initial force required to overcome static friction and other resistive forces in many electric motors
but still allows it to continue operating down to a force lower than that initial requirement.

## Mathematical formulation

The cart-pole *n*-link system has configuration $q = (x, \theta_1, \dots, \theta_n)$
with $\theta_i = 0$ hanging straight down and $\theta_i = \pi$ upright. The
equations of motion follow from the Euler–Lagrange equations and take the
standard manipulator form:

$$
M(q)\mathbf{\ddot{q}} + \mathbf{h}(q, \dot{q}) = B \mathbf{u}
$$

where $M$ is the configuration-dependent inertia matrix, $\mathbf{h}$ contains the velocity cross terms and the gravitational terms, and $B = [1, 0, \dots, 0]^\top$
since the cart is the only actuated coordinate.

The linearization is done by a finite difference of the state around a unstable equilibrium. The initial plan was to derive the 

The LQR cost is

$$
J = \int_0^\infty \left( z^\top Q z + u^\top R u \right) dt
$$

with [VERIFY: $Q = \mathrm{diag}(\ldots)$, $R = \ldots$]. Solving the algebraic
Riccati equation gives the optimal gain $K$ which would give $u = -Kx$.

## Repository layout

```
src/
  PendulumConfig.py  Pendulum class with class methods for configuring the input parameters and ICs
  PendulumSystem.py  Pendulum class with methods for calculating inertia, coriolis, gravity terms, finding the rhs of the diff eqn and then running the simulation
  controllers.py     LQR synthesis and deadband wrapper
  visualization.py   Manim animation tools
notebooks/
  01_single_pendulum.ipynb
  02_double_pendulum.ipynb
  03_triple_pendulum.ipynb
tests/
  test_energy_conservation.py
```

## Quick start

```bash
git clone https://github.com/andrewlastname/cartpole-lqr-control.git
cd cartpole-lqr-control
pip install -r requirements.txt

# Run the triple pendulum demo
python -m src.simulation --system triple

# Or open a notebook
jupyter notebook notebooks/03_triple_pendulum.ipynb
```

## What's next

- **Map regions of stability**: Use monte carlo style simulations to map out what the
  shape of the region in which stability is recoverable from the starting initial conditions.
  Expect this to have a elliptoid shape.
- **Reinforement Learning control**: use reinforcement learning to progressively
  train a model to control a single then double then triple pendulum. 
- **Swing-up controller**: energy-shaping or hybrid LQR to bring the system to
  upright from rest, then hand off to the stabilizing controller.
- **MPC comparison**: receding-horizon control on the same plant, with input
  constraints, to compare against LQR on transient performance and constraint
  handling.
- **Observer design**: Kalman filter for noisy measurements, leading to LQG.
- **Robustness analysis**: parameter uncertainty in masses and lengths,
  evaluated against the closed-loop spectrum.

## References

Russ Tedrake. Underactuated Robotics: Algorithms for Walking, Running, Swimming, Flying, and Manipulation (Course Notes for MIT 6.832). 
Downloaded on 05/20 from https://underactuated.csail.mit.edu/

Brunton, S. L. "Control Bootcamp." YouTube, University of Washington.
https://www.youtube.com/playlist?list=PLMrJAkhIeNNR20Mz-VpzgfQs5zrYi085m

## Acknowledgments

Manim animation scaffolding and debugging assistance provided by Claude
(Anthropic). The control design, dynamics derivation, and validation are
my own work.
