from abc import ABC, abstractmethod
from scipy.linalg import solve_continuous_are
import numpy as np

class Controller(ABC):
    """Base class. Every controller takes (t, y) and returns control input u."""
    @abstractmethod
    def __call__(self, t, y):
        """Return control input u for time t and state y."""
        pass

    def reset(self):
        """Override if your controller has internal state (integrators, observers)."""
        pass
    
    
class ZeroController(Controller):
    def __call__(self, t, y):
        return 0.0


class LQRController(Controller):
    def __init__(self, K, y_ref):
        self.K = np.atleast_2d(K)          # ensure 2D even for single inputs
        self.y_ref = np.asarray(y_ref)

    def __call__(self, t, y):
        u = -self.K @ (y - self.y_ref)
        return float(u[0]) if u.size == 1 else u

    @classmethod
    def design(cls, system, y_eq, Q, R, u_eq=0.0):
        """Compute K from system linearization and return a ready controller."""
        A, B = system.linearize(y_eq, u_eq)
        P = solve_continuous_are(A, B, Q, R)
        K = np.linalg.inv(R) @ B.T @ P
        return cls(K, y_eq)
    
class DeadbandActuator(Controller):
    def __init__(self, inner_controller, deadband_on, deadband_off):
        self.inner_controller = inner_controller
        self.deadband_on = deadband_on
        self.deadband_off = deadband_off
        self._controller_active = False

    def __call__(self, t, y):
        u = self.inner_controller(t, y)
        if abs(u) > self.deadband_off and abs(u) < self.deadband_on:
            self._controller_active = True
            return u*(abs(u) - self.deadband_off) / (self.deadband_on - self.deadband_off) # Scale u to be continuous at the deadband boundaries. Gives a linear decrease down to zero at the off point.
        elif abs(u) > self.deadband_on:
            self._controller_active = True
            return u
        else:
            self._controller_active = False
            return 0.0
