from dataclasses import dataclass, field
import numpy as np

@dataclass
class PendulumConfig:
  n: int
  lengths: np.ndarray
  masses: np.ndarray
  is_cart: bool = False
  gravity: float = 9.81
  
  def __post_init__(self):
    self.lengths = np.asarray(self.lengths, dtype=float)
    self.masses = np.asarray(self.masses, dtype=float)
    assert len(self.lengths) == self.n, (
      f"lengths must have length n={self.n}, got {len(self.lengths)}")
    assert len(self.masses) in (self.n, self.n + 1), (
      f"masses must have length n = {self.n} or {self.n+1}, got {len(self.masses)}")
      
  @classmethod
  def default_params(cls, n, is_cart: bool = False):
    lengths = np.ones(n)
    masses = np.concatenate([[10], np.ones(n)]) if is_cart else np.ones(n)
    return cls(n=n, lengths=lengths, masses=masses, is_cart=is_cart)
  
  def small_angle_ICs(self, eps: float = .1):
    """Small initial angles for a small angle approximation sim measured from the hanging vertical position"""
    theta_init = eps*np.ones(self.n)
    theta_dot_init = np.zeros(self.n)
    return np.concatenate([[0], theta_init, [0], theta_dot_init]) if self.is_cart else np.concatenate([theta_init, theta_dot_init])
  
  def chaotic_ICs(self, eps: float = .01):
    """Large angles for a chaotic sim measured from the hanging vertical position"""
    theta_init = (np.pi-eps)*np.ones(self.n)
    theta_dot_init = np.zeros(self.n)
    return np.concatenate([[0], theta_init, [0], theta_dot_init]) if self.is_cart else np.concatenate([theta_init, theta_dot_init])
  
  def random_ICs(self, eps: float = .01):
    """Returns a random set of initial conditions for the pendulum sim measured from the hanging vertical position"""
    theta_init = np.random.uniform(-np.pi, np.pi, size=self.n)
    theta_dot_init = np.random.uniform(-1, 1, size=self.n)
    x_init = [np.random.uniform(-1, 1)] if self.is_cart else None
    x_dot_init = [0] if self.is_cart else None #Keeping zero for the purpose of preventing the cart from flying off in the random sims, but this could be changed to be random as well.
    return np.concatenate([x_init, theta_init, x_dot_init, theta_dot_init]) if self.is_cart else np.concatenate([theta_init, theta_dot_init])