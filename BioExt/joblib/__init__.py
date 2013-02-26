
try:
    from joblib import Parallel, delayed, __version__
    assert __version__.find('lazypar') > 0
except (AssertionError, ImportError):
    import warnings

    warnings.warn(
        'joblib not present or incompatible, falling back to a single-threaded implementation'
        )

    def delayed(function):
        return function

    class Parallel:
        def __init__(self, n_jobs=None, verbose=None, pre_dispatch=None):
            pass

        def lazy(self, iterable):
            return iterable

        def __call__(self, iterable):
            return iterable
