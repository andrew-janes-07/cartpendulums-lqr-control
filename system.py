import numpy as np
from scipy.integrate import solve_ivp
from src.config import PendulumConfig
class PendulumSystem:
  def __init__(self, n, is_cart=False,cfg = None, ICs=None):
    self.n = n
    self.is_cart = is_cart
    self.cfg = cfg if cfg is not None else PendulumConfig.default_params(n, is_cart=self.is_cart)
    self.init_state = ICs if ICs is not None else self.cfg.small_angle_ICs()
    self.t = None
    self.state = None
    self.theta_init = self.init_state[1:self.n+1] if self.is_cart else self.init_state[:self.n] # Initial Angles
    self.theta_dot_init = self.init_state[self.n+2:] if self.is_cart else self.init_state[self.n:] # Initial Angular Velocities
    
    self.m_T = self.cfg.masses.sum() # Only used for the cart pendulum system, gives the system mass either way though
    if self.n == 1:
      self.mu = np.array([self.cfg.masses[1:].sum()]) if self.is_cart else np.array([self.cfg.masses[0:].sum()])
    elif self.n == 2:
      self.mu = np.array([self.cfg.masses[1:].sum(), self.cfg.masses[2:].sum()]) if self.is_cart else np.array([self.cfg.masses[0:].sum(), self.cfg.masses[1:].sum()])
    else: 
      self.mu = np.array([self.cfg.masses[1:].sum(), self.cfg.masses[2:].sum(), self.cfg.masses[3:].sum()]) if self.is_cart else np.array([self.cfg.masses[0:].sum(), self.cfg.masses[1:].sum(), self.cfg.masses[2:].sum()])
      
    self.control_history = None # This will be used to store the control history for energy calculations when a controller is used in the simulate function. It will be initialized as an empty list when a controller is used in the simulate function and then appended to at each time step within the solve_eoms function.
  
  @property
  def theta(self):
      if self.state is None:
        return self.theta_init
      else:
        return self.state[1:self.n+1] if self.is_cart else self.state[:self.n]
      
      
  @property
  def theta_dot(self):
      if self.state is None:
        return self.theta_dot_init
      else:
        return self.state[self.n+2:] if self.is_cart else self.state[self.n:]
    
    
  def inertia_matrix(self, theta = None, theta_dot = None):
    """ 
  Computes the inertia matrix for the pendulum equations of motion.
    
  Parameters:
  theta : array, optional
      Angles off the bottom vertical. Takes inputs for testing purposes or calling separately from the simulator. If None, will use self.theta. Otherwise it defaults to the current state of the pendulum system.
    
  Returns:
  M : array
      Inertia Matrix
  """
    if theta is None: theta = self.theta
    n = self.n
    L = self.cfg.lengths
    mu = self.mu
      
    if self.is_cart:
      # Calculates the inertia matrix for n pendulums on a cart
      M = np.zeros((self.n+1, self.n+1))
      for i in range(self.n + 1):
        for j in range(self.n + 1):
          if i == 0 and j == 0:
            M[i, j] = self.m_T
          elif i == j:
            M[i, j] = mu[i-1]*L[i-1]**2
          elif i == 0:
            M[i, j] = mu[j-1]*L[j-1]*np.cos(theta[j-1])
          elif j > i:
            M[i, j] = mu[j-1]*L[i-1]*L[j-1]*np.cos(theta[i-1]-theta[j-1])
          else:
            M[i, j] = M[j, i] # Symmetry terms
      return M
        
    else:
      # Calculates the inertia matrix for n pendulums without a cart
      M = np.zeros((self.n, self.n))
      for i in range(self.n):
        for j in range(self.n):
          if i == j:
            M[i, j] = mu[i]*L[i]**2
          else:
            M[i, j] = mu[max(i, j)]*L[i]*L[j]*np.cos(theta[i]-theta[j])
      return M
  
  
  def nonlinear_terms_vector(self, theta = None, theta_dot = None):
    """ 
    Computes the coriolis and gravity terms for the pendulum equations of motion.
    
    Parameters:
    theta : array, optional
        Angles off the bottom vertical. Takes inputs for testing purposes or calling separately from the simulator. If None, will use self.theta. Otherwise it defaults to the current state of the pendulum system.
    theta_dot : array, optional
        Angular velocities in radians per second. Takes inputs for testing purposes or calling separately from the simulator. If None, will use self.theta. Otherwise it defaults to the current state of the pendulum system. 
    
    Returns:
    C + G : array
        Nonlinear terms vector.
    """
    
    if theta is None: theta = self.theta
    if theta_dot is None: theta_dot = self.theta_dot
    n = self.n
    L = self.cfg.lengths
    mu = self.mu
    g = self.cfg.gravity
    
    if self.is_cart:
      # Calculates the nonlinear terms for n pendulums on a cart
      C = np.zeros(self.n+1)
      G = np.zeros(self.n+1)
      C[0] = -np.sum(mu*L*np.sin(theta)*theta_dot**2)
      G[0] = 0
      offset = 1
    else: 
      C = np.zeros(self.n)
      G = np.zeros(self.n)
      offset = 0
    for i in range(self.n):
      G[i + offset] = mu[i]*L[i]*g*np.sin(theta[i])
      for j in range(0, self.n):
        if j != i:
          C[i+offset] += mu[max(i,j)]*L[i]*L[j]*np.sin(theta[i]-theta[j])*theta_dot[j]**2
    return C+G
  
  
  def solve_eoms(self, u):
    """ 
    Calls inertia matrix and nonlinear term functions and then solves for q_ddot at that specific state
    """
    M = self.inertia_matrix()
    h = self.nonlinear_terms_vector()
    if self.is_cart:
      h[0] += u
    q_ddot = -np.linalg.solve(M, h)
    return q_ddot
  
  
  def _rhs(self, t, y, u=0.0):
    """State derivative function for solve_ivp with control implemented as a force applied to the cart if is_cart is True. This function is not meant to be called directly, but rather is used within the simulate function of the class. It takes in the current state and time and outputs the state derivatives for the next step in the simulation."""
    self.state = y
    q_ddot = self.solve_eoms(u)
    if self.is_cart:
      # y = [x, theta_1, theta_2, ..., x_dot, theta_1_dot, theta_2_dot, ...]
      # y_dot = [x_dot, theta_1_dot, theta_2_dot, ..., x_ddot, theta_1_ddot, theta_2_ddot, ...]
      q_dot = y[self.n+1:]
      return np.concatenate([q_dot, q_ddot])
    else:
      # y = [theta_1, theta_2, ..., theta_1_dot, theta_2_dot, ...]
      # y_dot = [theta_1_dot, theta_2_dot, ..., theta_1_ddot, theta_2_ddot, ...]
      q_dot = y[self.n:]
      return np.concatenate([q_dot, q_ddot])
    
    
  def simulate(self, solver = 'DOP853', t_span = (0, 10), dt = .01, controller=None):
    """ 
    Simulates the pendulum system using DOP853 (default) and returns the time points and state history.
    
    Parameters:
    solver : str
        The solver to use for the simulation.
    t_span : tuple
        Start and end times for the simulation.
    dt : float
        Time step for evaluation points.
    
    Returns:
    sim.t : array
        Time points of the simulation.
    sim.y : array
        State history of the simulation, with shape (n_states, n_time_points).
    """
    if controller is not None:
        controller.reset()           # clear any internal state from prior runs
    u_fn = controller if controller is not None else (lambda t, y: 0.0)
    
    eval_points = int((t_span[1] - t_span[0]) / dt) + 1
    t_eval = np.linspace(t_span[0], t_span[1], eval_points)
    sim = solve_ivp(lambda t, y: self._rhs(t, y, u_fn(t, y)), t_span=t_span, y0=self.init_state, t_eval=t_eval, method=solver, rtol=1e-10, atol=1e-12)
    self.sim_t = sim.t # Creates an accessible t vector for energy calculations within the class 
    self.sim_y = sim.y # Creates an accessible state history vector for energy calculations within the class
    
    # Recompute control input on the clean output grid
    self.control_history = np.array([
        u_fn(t, y) for t, y in zip(sim.t, sim.y.T)
    ])
    
    return sim
    
  def linearize(self, y_eq, u_eq=0.0, eps=1e-6):
    """u_eq gives default control input for linearization, which is zero for the pendulum on a cart system. y_eq gives the state for linearization."""
    
    saved_state = self.state    # save
    try:
        n_y = len(y_eq)
        A = np.zeros((n_y, n_y))
        for i in range(n_y):
            y_plus = y_eq.copy();  y_plus[i] += eps
            y_minus = y_eq.copy(); y_minus[i] -= eps
            A[:, i] = (self._rhs(0.0, y_plus, u_eq) - self._rhs(0.0, y_minus, u_eq)) / (2 * eps)
        
        if self.is_cart:
            f_plus = self._rhs(0.0, y_eq, u_eq + eps)
            f_minus = self._rhs(0.0, y_eq, u_eq - eps)
            B = ((f_plus - f_minus) / (2 * eps))[:, None]
        else:
            B = np.zeros((n_y, 0))
        
        return A, B
    finally:
        self.state = saved_state    # restore
    
  def energy_history(self):
    """ 
    Calculate the energy of the system at each time point in the simulation using E = T + V. This function should be called after simulate to access the time points and state history for energy calculations.
    
    Returns:
    E : array
        Energy of the system at each time point in the simulation.
    """
    if self.sim_t is None or self.sim_y is None:
      raise ValueError("No simulation data found. Please run simulate() before calling state_energy().")
    n = self.n
    n_t = self.sim_t.shape[0]
    g = self.cfg.gravity
    L = self.cfg.lengths
    mu = self.mu
    
    if self.is_cart:
      theta_history = self.sim_y[1:n+1]
      q_dot_history = self.sim_y[n+1:]
      dim = n+1
    else:
      theta_history = self.sim_y[:n]
      q_dot_history = self.sim_y[n:]
      dim = n
    
    U = -g*(mu*L)@np.cos(theta_history) # .reshape([n_t])
    M_history = np.zeros((n_t, dim, dim))
    for i in range(n_t):
      M_history[i, :, :] = self.inertia_matrix(theta=theta_history[:, i])
    T = .5*np.einsum('it,tij,jt->t', q_dot_history, M_history, q_dot_history, optimize=True)
    return T.T, U.T, (T+U).T