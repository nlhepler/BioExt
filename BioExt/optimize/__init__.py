"""
This was all lifted from scipy master on June 6, 2013.
The license for this code is therefore different from everything else,
and is the standard 3-clause BSD license found in LICENSE.txt in the
same directory as this file.
"""
from __future__ import division, print_function

import numpy

from BioExt.optimize._linesearch import line_search_wolfe2


__all__ = ['minimize_bfgs']


_epsilon = numpy.sqrt(numpy.finfo(float).eps)
# standard status messages of optimizers
_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.'}


class Result(dict):
    """ Represents the optimization result.

    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess : ndarray
        Values of objective function, Jacobian and Hessian (if available).
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in self.items()])
        else:
            return self.__class__.__name__ + "()"


def wrap_function(function, args):
    ncalls = [0]
    if function is None:
        return ncalls, None

    def function_wrapper(*wrapper_args):
        ncalls[0] += 1
        return function(*(wrapper_args + args))

    return ncalls, function_wrapper


def vecnorm(x, ord=2):
    if ord == numpy.Inf:
        return numpy.amax(numpy.abs(x))
    elif ord == -numpy.Inf:
        return numpy.amin(numpy.abs(x))
    else:
        return numpy.sum(numpy.abs(x) ** ord, axis=0) ** (1.0 / ord)


def approx_fprime(xk, f, epsilon, *args):
    """Finite-difference approximation of the gradient of a scalar function.

    Parameters
    ----------
    xk : array_like
        The coordinate vector at which to determine the gradient of `f`.
    f : callable
        The function of which to determine the gradient (partial derivatives).
        Should take `xk` as first argument, other arguments to `f` can be
        supplied in ``*args``.  Should return a scalar, the value of the
        function at `xk`.
    epsilon : array_like
        Increment to `xk` to use for determining the function gradient.
        If a scalar, uses the same finite difference delta for all partial
        derivatives.  If an array, should contain one value per element of
        `xk`.
    \*args : args, optional
        Any other arguments that are to be passed to `f`.

    Returns
    -------
    grad : ndarray
        The partial derivatives of `f` to `xk`.

    See Also
    --------
    check_grad : Check correctness of gradient function against approx_fprime.

    Notes
    -----
    The function gradient is determined by the forward finite difference
    formula::

                 f(xk[i] + epsilon[i]) - f(xk[i])
        f'[i] = ---------------------------------
                            epsilon[i]

    The main use of `approx_fprime` is in scalar function optimizers like
    `fmin_bfgs`, to determine numerically the Jacobian of a function.

    Examples
    --------
    >>> from scipy import optimize
    >>> def func(x, c0, c1):
    ...     "Coordinate vector `x` should be an array of size two."
    ...     return c0 * x[0]**2 + c1*x[1]**2

    >>> x = np.ones(2)
    >>> c0, c1 = (1, 200)
    >>> eps = np.sqrt(np.finfo(np.float).eps)
    >>> optimize.approx_fprime(x, func, [eps, np.sqrt(200) * eps], c0, c1)
    array([   2.        ,  400.00004198])

    """
    f0 = f(*((xk,) + args))
    grad = numpy.zeros((len(xk),), float)
    ei = numpy.zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei
        grad[k] = (f(*((xk + d,) + args)) - f0) / d[k]
        ei[k] = 0.0

    return grad


class _LineSearchError(RuntimeError):
    pass


def _line_search_wolfe2(f, fprime, xk, pk, gfk, old_fval, old_old_fval,
                        **kwargs):
    """
    Same as line_search_wolfe1, but fall back to line_search_wolfe2 if
    suitable step length is not found, and raise an exception if a
    suitable step length is not found.

    Raises
    ------
    _LineSearchError
        If no suitable step size is found

    """
    ret = line_search_wolfe2(f, fprime, xk, pk, gfk,
                             old_fval, old_old_fval,
                             **kwargs)

    if ret[0] is None:
        raise _LineSearchError()

    return ret


def minimize_bfgs(fun, x0, args=(), jac=None, callback=None,
                  tol=1e-5, norm=numpy.Inf, eps=_epsilon, maxiter=None,
                  disp=False, return_all=False):
    """
    Minimization of scalar function of one or more variables using the
    BFGS algorithm.

    Options for the BFGS algorithm are:
        disp : bool
            Set to True to print convergence messages.
        maxiter : int
            Maximum number of iterations to perform.
        tol : float
            Gradient norm must be less than `tol` before successful
            termination.
        norm : float
            Order of norm (Inf is max, -Inf is min).
        eps : float or ndarray
            If `jac` is approximated, use this value for the step size.
    """
    f = fun
    fprime = jac
    epsilon = eps
    retall = return_all

    x0 = numpy.asarray(x0).flatten()
    if x0.ndim == 0:
        x0.shape = (1,)
    if maxiter is None:
        maxiter = len(x0) * 200
    func_calls, f = wrap_function(f, args)
    if fprime is None:
        grad_calls, myfprime = wrap_function(approx_fprime, (f, epsilon))
    else:
        grad_calls, myfprime = wrap_function(fprime, args)
    gfk = myfprime(x0)
    k = 0
    N = len(x0)
    I = numpy.eye(N, dtype=int)
    Hk = I
    old_fval = f(x0)
    old_old_fval = None
    xk = x0
    if retall:
        allvecs = [x0]
    sk = [2 * tol]
    warnflag = 0
    gnorm = vecnorm(gfk, ord=norm)
    while (gnorm > tol) and (k < maxiter):
        pk = -numpy.dot(Hk, gfk)
        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                _line_search_wolfe2(f, myfprime, xk, pk, gfk,
                                    old_fval, old_old_fval)
        except _LineSearchError:
            # Line search failed to find a better solution.
            warnflag = 2
            break

        xkp1 = xk + alpha_k * pk
        if retall:
            allvecs.append(xkp1)
        sk = xkp1 - xk
        xk = xkp1
        if gfkp1 is None:
            gfkp1 = myfprime(xkp1)

        yk = gfkp1 - gfk
        gfk = gfkp1
        if callback is not None:
            callback(xk)
        k += 1
        gnorm = vecnorm(gfk, ord=norm)
        if (gnorm <= tol):
            break

        if not numpy.isfinite(old_fval):
            # We correctly found +-Inf as optimal value, or something went
            # wrong.
            warnflag = 2
            break

        try:  # this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (numpy.dot(yk, sk))
        except ZeroDivisionError:
            rhok = 1000.0
            if disp:
                print("Divide-by-zero encountered: rhok assumed large")
        if numpy.isinf(rhok):  # this is patch for numpy
            rhok = 1000.0
            if disp:
                print("Divide-by-zero encountered: rhok assumed large")
        A1 = I - sk[:, numpy.newaxis] * yk[numpy.newaxis, :] * rhok
        A2 = I - yk[:, numpy.newaxis] * sk[numpy.newaxis, :] * rhok
        Hk = numpy.dot(A1, numpy.dot(Hk, A2)) + (rhok * sk[:, numpy.newaxis] *
                                                 sk[numpy.newaxis, :])

    fval = old_fval
    if warnflag == 2:
        msg = _status_message['pr_loss']
        if disp:
            print("Warning: " + msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])

    elif k >= maxiter:
        warnflag = 1
        msg = _status_message['maxiter']
        if disp:
            print("Warning: " + msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])
    else:
        msg = _status_message['success']
        if disp:
            print(msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])

    result = Result(fun=fval, jac=gfk, hess=Hk, nfev=func_calls[0],
                    njev=grad_calls[0], status=warnflag,
                    success=(warnflag == 0), message=msg, x=xk)
    if retall:
        result['allvecs'] = allvecs
    return result
