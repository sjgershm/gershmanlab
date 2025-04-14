from scipy.stats import ttest_1samp, ttest_ind
from scipy.integrate import quad
import numpy as np

def bayesian_ttest(x, y=None, r=1):
    """
    Bayes Factor t-test (one-sample or two-sample) following Rouder et al. (2009)

    Parameters:
    - x: array-like, sample data or group 1
    - y: array-like or None. If None, perform one-sample test; else, two-sample test with group 2
    - r: float, Cauchy scale factor (default = 1)

    Returns:
    - B: Bayes factor favoring the null hypothesis
    """
    x = np.asarray(x)
    if y is None:
        # One-sample t-test
        N = len(x)
        v = N - 1
        t_stat, _ = ttest_1samp(x, 0)
    else:
        y = np.asarray(y)
        N1 = len(x)
        N2 = len(y)
        N = (N1 * N2) / (N1 + N2)  # effective sample size for independent samples
        v = N1 + N2 - 2
        t_stat, _ = ttest_ind(x, y, equal_var=True)

    t = t_stat
    numerator = (1 + (t**2) / v) ** (-(v + 1) / 2)

    def integrand(g):
        denom = (1 + N * g * r**2)
        term1 = denom**(-0.5)
        term2 = (1 + (t**2) / (denom * v)) ** (-(v + 1) / 2)
        term3 = (2 * np.pi)**(-0.5) * g**(-1.5) * np.exp(-1 / (2 * g))
        return term1 * term2 * term3

    denominator, _ = quad(integrand, 0, np.inf)
    B = numerator / denominator
    return B
