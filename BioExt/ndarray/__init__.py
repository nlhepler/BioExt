
import numpy as np


__all__ = [
    'symndarray',
    'symzeros'
    ]


class symndarray(np.ndarray):

    def __setitem__(self, indices, value):
        try:
            i, j = indices[:2]
            indices_ = (j, i) + indices[2:]
            np.ndarray.__setitem__(self, indices_, value)
            np.ndarray.__setitem__(self, indices, value)
        except:
            self[indices, :] = value


def symzeros(*args, **kwargs):
    return np.zeros(*args, **kwargs).view(symndarray)
