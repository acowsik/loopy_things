# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import random
import pickle
from typing import List, Tuple
import numpy as np
cimport numpy as np

cdef class Ising:
    cdef public unsigned char[:] state
    cdef public long x_size
    cdef public long y_size
    cdef public long N
    cdef long long rand_index
    cdef double[:] random_nums
    cdef double temperature

    def __init__(self, log_x_size, log_y_size, temperature):
        # initialize a state of size which is a power of 2 for future bit twiddling hacks
        self.x_size = 2**log_x_size
        self.y_size = 2**log_y_size
        self.N = self.x_size * self.y_size

        self.state = np.random.randint(2, size=self.N, dtype=np.uint8)
        self.temperature = temperature

        # generate a bunch of random numbers initially, I have found this is more 
        # efficient somehow
        self.rand_index = 0
        self.rand_size = 16 * self.N
        self.random_spins = np.random.random(size=self.rand_size)
        

    def step(self):
        cdef int spin = self.rand_index
        self.rand_index += 1

        # refresh the random numbers if needed
        # works because rand_size is a power of 2
        if self.rand_index & self.rand_size: 
            self.random_spins = np.random.random(size=self.rand_size)
            self.rand_index = 0
        


