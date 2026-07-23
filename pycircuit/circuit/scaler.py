from abc import ABC, abstractmethod

class Scaler(ABC):
    """
    Abstract Strategy Interface for Jacobian Matrix Scaling/Equilibration.
    Scaling improves the numerical condition number of the Jacobian, preventing
    precision loss during the linear solve step of the Newton-Raphson iteration.
    """
    @abstractmethod
    def scale(self, J, F, toolkit):
        """
        Scales the Jacobian matrix J and the right-hand-side vector F.
        Returns:
            J_scaled, F_scaled, scale_vector
        """
        pass
        
    @abstractmethod
    def unscale_solution(self, dx, scale_vector, toolkit):
        """
        Some scaling methods (like column scaling) require unscaling the 
        resulting solution vector (dx). Row scaling usually does not.
        """
        pass


class NoneScaler(Scaler):
    """Dummy scaler that performs no scaling."""
    def scale(self, J, F, toolkit):
        return J, F, None
        
    def unscale_solution(self, dx, scale_vector, toolkit):
        return dx


class RowMaxScaler(Scaler):
    """
    Divides each row of J and F by the maximum absolute value in that row.
    Extremely fast, prevents massive values from dominating the solve.
    """
    def scale(self, J, F, toolkit):
        # Calculate maximum absolute value per row
        row_max = toolkit.max(abs(J), axis=1)
        
        # Prevent division by zero for empty rows (or rows with 0 max)
        row_max = toolkit.where(row_max < 1e-20, 1.0, row_max)
        
        # Reshape for broadcasting
        # This divides each element in row i by row_max[i]
        J_scaled = J / row_max[:, None]
        F_scaled = F / row_max
        
        return J_scaled, F_scaled, None # No solution unscaling needed for row scaling
        
    def unscale_solution(self, dx, scale_vector, toolkit):
        return dx


class RowL2Scaler(Scaler):
    """
    Divides each row by its L2-Norm (Euclidean length).
    Provides smoother statistical balancing than Max scaling.
    """
    def scale(self, J, F, toolkit):
        # Calculate L2 norm per row: sqrt(sum(x^2))
        row_l2 = toolkit.sqrt(toolkit.sum(J**2, axis=1))
        
        # Prevent division by zero
        row_l2 = toolkit.where(row_l2 < 1e-20, 1.0, row_l2)
        
        J_scaled = J / row_l2[:, None]
        F_scaled = F / row_l2
        
        return J_scaled, F_scaled, None
        
    def unscale_solution(self, dx, scale_vector, toolkit):
        return dx


class SinkhornKnoppScaler(Scaler):
    """
    Iterative equilibration that scales both rows and columns to approach a 
    doubly stochastic matrix (all row and column sums equal 1).
    Minimizes condition number robustly for poorly conditioned circuits.
    """
    def __init__(self, max_iter=5):
        self.max_iter = max_iter

    def scale(self, J, F, toolkit):
        J_s = J * 1.0
        F_s = F * 1.0
        
        # We must track the accumulated column scaling to unscale dx later
        c_scale_total = toolkit.ones(J.shape[1])
        
        for _ in range(self.max_iter):
            # 1. Row scaling
            r_sums = toolkit.sum(abs(J_s), axis=1)
            r_sums = toolkit.where(r_sums < 1e-20, 1.0, r_sums)
            J_s = J_s / r_sums[:, None]
            F_s = F_s / r_sums
            
            # 2. Column scaling
            c_sums = toolkit.sum(abs(J_s), axis=0)
            c_sums = toolkit.where(c_sums < 1e-20, 1.0, c_sums)
            J_s = J_s / c_sums[None, :]
            
            # Accumulate the column divisors
            c_scale_total = c_scale_total * c_sums
            
        return J_s, F_s, c_scale_total
        
    def unscale_solution(self, dx, scale_vector, toolkit):
        # We divided the columns of J by scale_vector, which means the variables
        # were implicitly multiplied by scale_vector. To recover the true dx, 
        # we must divide the solved dx by the scale_vector.
        return dx / scale_vector
