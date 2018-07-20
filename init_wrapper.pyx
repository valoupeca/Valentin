# -*- coding: utf-8 -*-
"""
Cython linker with C solver
"""

# Author: Remi Flamary <remi.flamary@unice.fr>
#
# License: MIT License

import numpy as np
cimport numpy as np

cimport cython

import warnings


cdef extern from "Dynamic_approach.h":
 mat Dynamic_calcul(vector<std::string> landmark)



@cython.boundscheck(False)
@cython.wraparound(False)
def Dynamic_calcul(np.array[string, ndim = 1, mode="c"]):

cdef np.array[string, ndim = 1, mode="c"]) landmark = np.zeros(20);



# calling the function
cdef mat result_code = Dynamic_calcul(landmark)

return result_code