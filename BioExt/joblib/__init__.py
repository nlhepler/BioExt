
try:
    try:
        from joblib import Parallel, delayed
    except ImportError:
        from sklearn.externals.joblib import Parallel, delayed
except ImportError:
    def delayed(function):
        return function
    class Parallel:
        def __init__(self, n_jobs=None, verbose=None, pre_dispatch=None):
            pass
        def __call__(self, iterable):
            return iterable
