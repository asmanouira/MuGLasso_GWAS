import spams 
from sklearn import base
from sklearn.utils import validation
from sklearn.utils import estimator_checks
from pre_data import input_matrix
from gsroptim.sgl_fast import bcd_fast
from gsroptim.sgl_tools import build_lambdas, precompute_norm

class MySpamsGroupLasso(base.BaseEstimator, base.RegressorMixin):
    """
    This class implements a group lasso that is fitted using SPAMS.
    """
    def __init__(self, lambda1=1., groups=None):
        self.lambda1 = lambda1
        self.groups = groups
                 
    def fit(self, X, y):     
        # Check that X and y have correct shape
        X, y = validation.check_X_y(X, y)
        
        self.X_ = X
        self.y_ = y     
                
        # Fitting 
        #X_asfa = np.asfortranarray(self.X_)
        X_sparse = sparse.csc_matrix(self.X_)        
        y_asfa = np.asfortranarray(self.y_.reshape((self.y_.shape[0], 1)))
        w_init = np.zeros((X_sparse.shape[1], 1), order="F")
        w = spams.fistaFlat(y_asfa, X_sparse, W0=w_init, loss='logistic', regul='group-lasso-l2', groups=self.groups, 
                            lambda1=self.lambda1)
        self.coef_ = w.reshape((w.shape[0], ))

        # Return the classifier
        return self
    
    def predict(self, X):
        # Check that fit has been called
        validation.check_is_fitted(self) 
        
        # Input validation
        X = validation.check_array(X)
        
        # Prediction
        return np.dot(X, self.coef_)

    def gap_safe_screening_rules (X, y, size_groups, omega, lambdas=None, tau=0, lambda2=0,
             beta_init=None, screen='GAPSAFE'):


        n_groups = len(size_groups)
        g_start = np.cumsum(size_groups, dtype=np.intc) - size_groups[0]

        if lambdas is None:
            lambdas, _ = build_lambdas(X, y, omega, size_groups, g_start)

        # Useful precomputation
        norm_X, norm_X_g, nrm2_y = precompute_norm(X, y, size_groups, g_start)
        tol = eps * nrm2_y  # duality gap tolerance

        n_lambdas = len(lambdas)
        n_samples, n_features = X.shape
        lambda_max = lambdas[0]

        # Fortran-contiguous array are used to avoid useless copy of the data.
        X = np.asfortranarray(X)
        y = np.asfortranarray(y)
        size_groups = np.asfortranarray(size_groups, dtype=np.intc)
        norm2_X = np.asfortranarray(norm_X ** 2)
        norm2_X_g = np.asfortranarray(norm_X_g ** 2)
        omega = np.asfortranarray(omega)
        g_start = np.asfortranarray(g_start, dtype=np.intc)

        if beta_init is None:
            beta_init = np.zeros(n_features, order='F')
        else:
            beta_init = np.asfortranarray(beta_init)

        betas = np.zeros((n_features, n_lambdas), order='F')
        residual = np.asfortranarray(y - np.dot(X, beta_init))
        XTR = np.asfortranarray(np.dot(X.T, residual))
        dual_scale = lambda_max  # good iif beta_init = 0

        gaps = np.ones(n_lambdas)
        screening_sizes_features = np.zeros(n_lambdas)
        screening_sizes_groups = np.zeros(n_lambdas)
        n_iters = np.zeros(n_lambdas)

        disabled_features = np.zeros(n_features, dtype=np.intc, order='F')
        disabled_groups = np.zeros(n_groups, dtype=np.intc, order='F')

        active_ws = False
        strong_ws = False

        for t in range(n_lambdas):

            if t == 0:
                lambda_prec = lambda_max
            else:
                lambda_prec = lambdas[t - 1]

            if strong_active_warm_start:
                strong_ws = True

            if gap_active_warm_start:
                active_ws = (screening_sizes_features[t] < n_features or
                             screening_sizes_groups[t] < n_groups)

            if strong_ws or active_ws:

                bcd_fast(X, y, beta_init, XTR, residual, dual_scale, omega,
                         n_samples, n_features, n_groups, size_groups, g_start,
                         norm2_X, norm2_X_g, nrm2_y, tau, lambdas[t],
                         lambda_prec, lambda2, max_iter, f, tol,
                         screen, disabled_features, disabled_groups,
                         wstr_plus=active_ws, strong_warm_start=strong_ws)

            model = bcd_fast(X, y, beta_init, XTR, residual, dual_scale, omega,
                             n_samples, n_features, n_groups, size_groups, g_start,
                             norm2_X, norm2_X_g, nrm2_y, tau, lambdas[t],
                             lambda_prec, lambda2, max_iter, f, tol,
                             screen, disabled_features, disabled_groups,
                             wstr_plus=0, strong_warm_start=0)

            dual_scale, gaps[t], _, _, n_iters[t] = model

            betas[:, t] = beta_init.copy()

            if verbose and abs(gaps[t]) > tol[t]:
                warnings.warn('Solver did not converge after '
                              '%i iterations: dual gap: %.3e'
                              % (max_iter, gaps[t]), ConvergenceWarning)

        return (betas, gaps, n_iters, screening_sizes_groups,
                screening_sizes_features)

