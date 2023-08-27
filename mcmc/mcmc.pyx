# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import random
import pickle
from typing import List, Tuple
import numpy as np
cimport numpy as np
from math import exp

def is_power_of_2(n):
    while n > 1:
        if n & 1:
            return False
        
        n = n >> 1
    
    return True

    
cdef class RandomLongs:
    cdef public long size
    cdef long[:] array
    cdef long rand_index
    cdef public long maximum

    def __init__(self, long size, long maximum):
        assert is_power_of_2(size)
        self.maximum = maximum
        self.size = size
        self.rand_index = 0
        self.refresh()

    def refresh(self):
        self.array = np.random.randint(0, self.maximum, size=self.size)
    
    cpdef long next(self):
        cdef long out = self.array[self.rand_index]
        self.rand_index += 1

        # refresh the random numbers if needed
        # works because rand_size is a power of 2
        if self.rand_index & self.size: 
            self.refresh()
            self.rand_index = 0
        
        return out

cdef class RandomDoubles:
    cdef public long size
    cdef double[:] array
    cdef long rand_index
    def __init__(self, long size):
        assert is_power_of_2(size)
        self.size = size
        self.rand_index = 0
        self.refresh()

    def refresh(self):
        self.array = np.random.random(size=self.size)
    
    cpdef double next(self):
        cdef double out = self.array[self.rand_index]
        self.rand_index += 1

        # refresh the random numbers if needed
        # works because rand_size is a power of 2
        if self.rand_index & self.size: 
            self.refresh()
            self.rand_index = 0
        
        return out


cdef class Ising:
    cdef public char[:] state
    cdef public long x_size
    cdef public long y_size
    cdef public long N
    cdef unsigned long Nm1
    cdef RandomDoubles random_nums
    cdef RandomLongs random_spins
    cdef double beta
    cdef double [:]precomputed_probabilities
    cdef long long spin_marker_value
    cdef long long spin_marks[:]

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
        # efficient
        self.random_spins = RandomLongs(16 * self.N, self.N)
        self.random_nums = RandomDoubles(16 * self.N)

        # precompute the exponentials to make things faster
        # about 4.5 times!
        self.precomputed_probabilities = np.exp(-2 * self.beta * (np.arange(9) - 4))

        # initialize the marker for the spins so that we don't have to backtrack 
        # and unmark spins
        self.spin_marker_value = 1
        self.spin_marks = np.zeros(self.N, dtype=np.int64)
        

    cpdef void step(self):
        """
        does one potential spin flip 
        """
        cdef unsigned long spin = self.random_spins.next()
        cdef double rand_val = self.random_nums.next()
        
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

    cpdef void wolff_sweep(self):
        cdef long spin = self.random_spins.next()
        



    cpdef long magnetization(self):
        """
        calculates the sum of all the spins
        """
        cdef int M = 0
        cdef int i
        for i in range(self.N):
            M += self.state[i]
        
        return M


