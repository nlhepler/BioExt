from libcpp.pair cimport pair
from libcpp.vector cimport vector

cdef extern from "math.hpp" namespace "math":
    double prob_background(double, double, int, int)

cdef extern from "rateclass.hpp" namespace "rateclass":
    cdef cppclass rateclass_t:
        rateclass_t(vector[pair[int, int]]&, int) except +
        void learn(double&, double&, vector[pair[double, double]]&)

cpdef double p_bg(double lg_rate, double lg_invrate, int cov, int k):
    return prob_background(lg_rate, lg_invrate, cov, k)

cdef class RateClass:
    cdef rateclass_t * _rc
    def __cinit__(self, data, factor):
        self._rc = new rateclass_t(data, factor)
    def __dealloc__(self):
        del self._rc
    def __call__(self):
        cdef double lg_L = 0.0, aicc = 0.0
        cdef vector[pair[double, double]] params
        self._rc.learn(lg_L, aicc, params)
        return lg_L, aicc, params
