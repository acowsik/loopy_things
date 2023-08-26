# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import random
import pickle
from typing import List, Tuple
import numpy as np
cimport numpy as np
from math import exp

cdef class Ising:
    cdef public char[:] state
    cdef public long x_size
    cdef public long y_size
    cdef public long N
    cdef unsigned long Nm1
    cdef long long rand_index
    cdef long rand_size
    cdef double[:] random_nums
    cdef long[:] random_spins
    cdef double beta
    cdef double [:]precomputed_probabilities

    def __init__(self, log_x_size, log_y_size, temperature):
        """
        Initialize the state with x_size 2**log_x_size and y_size 2**log_y_size and at the given temperature 
        """
        # initialize a state of size which is a power of 2 for future bit twiddling hacks
        self.x_size = 2**log_x_size
        self.y_size = 2**log_y_size
        self.N = self.x_size * self.y_size
        self.Nm1 = (<unsigned long>self.N) - 1

        self.state = np.random.randint(2, size=self.N, dtype=np.uint8) * 2 - 1
        self.beta = 1.0 / temperature

        # generate a bunch of random numbers initially, I have found this is more 
        # efficient somehow
        self.rand_index = 0
        self.rand_size = 16 * self.N
        self.random_spins = np.random.randint(0, self.N, size=self.rand_size)
        self.random_nums = np.random.random(size=self.rand_size)

        # precompute the exponentials to make things faster
        # about 4.5 times!
        self.precomputed_probabilities = np.exp(-2 * self.beta * (np.arange(9) - 4))
        

    cpdef void step(self):
        """
        does one potential spin flip 
        """
        cdef unsigned long spin = self.random_spins[self.rand_index]
        cdef double rand_val = self.random_nums[self.rand_index]
        self.rand_index += 1

        # refresh the random numbers if needed
        # works because rand_size is a power of 2
        if self.rand_index & self.rand_size: 
            self.random_spins = np.random.randint(0, self.N, size=self.rand_size)
            self.random_nums = np.random.random(size=self.rand_size)
            self.rand_index = 0
        
        # find the net spin of the neighbors
        # the "& self.Nm1" is equivalent to going modulo self.N
        cdef int neighboring_values_total = 0
        neighboring_values_total += self.state[(spin + 1) & self.Nm1]
        neighboring_values_total += self.state[(spin - 1) & self.Nm1]
        neighboring_values_total += self.state[(spin + self.x_size) & self.Nm1]
        neighboring_values_total += self.state[(spin - self.x_size) & self.Nm1]

        cdef double one_probability = 1/(1 + self.precomputed_probabilities[neighboring_values_total + 4])
        self.state[spin] = <unsigned char> (rand_val < one_probability) * 2 - 1
    
    cpdef void sweep(self):
        """
        does N potential spin flips
        """
        cdef int i
        for i in range(self.N):
            self.step()


    cpdef long magnetization(self):
        """
        calculates the sum of all the spins
        """
        cdef int M = 0
        cdef int i
        for i in range(self.N):
            M += self.state[i]
        
        return M


