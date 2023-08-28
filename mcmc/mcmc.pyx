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

cdef class Stack:
    cdef public long long[:] stack
    cdef public unsigned long long stack_ptr
    cdef public long long stack_size
    def __init__(self, size):
        self.stack = np.zeros(size, dtype=np.int64)
        self.stack_size = size
        self.stack_ptr = 0
    
    cpdef void push(self, long long x):
        # doesn't check if we buffer overrun because that's done automatically
        self.stack[self.stack_ptr] = x
        self.stack_ptr += 1
    
    cpdef long long pop(self):
        # again does not check if we buffer overrun
        self.stack_ptr -= 1
        return self.stack[self.stack_ptr]

    
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

    cpdef refresh(self):
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
    cdef public unsigned long x_size
    cdef public unsigned long y_size
    cdef public long N
    cdef unsigned long Nm1
    cdef RandomDoubles random_nums
    cdef RandomLongs random_spins
    cdef double beta
    cdef double[:] precomputed_probabilities
    cdef long long spin_marker_value
    cdef long long[:] spin_marks
    cdef Stack wolff_stack

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

        # initialize the stack for the wolff algorithm so we don't
        # have to make a new one each time
        # the factor of 4 is because a spin might be added up to 4 times 
        # by each of its neighbors and 1 for safety
        self.wolff_stack = Stack(self.N * 4 + 1)
        

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

    cpdef long wolff_sweep(self):
        # from "Collective Monte Carlo Updating for Spin Systems" Ulli Wolff
        # choose a random lattice site
        cdef unsigned long spin = self.random_spins.next()
        self.wolff_stack.push(spin)
        self.spin_marker_value += 1
        cdef double probability = 1 - exp(-2 * self.beta)
        cdef long neighbor_spin

        cdef long spins_examined = 0

        while self.wolff_stack.stack_ptr > 0:
            spin = self.wolff_stack.pop()
            if self.spin_marks[spin] != self.spin_marker_value:
                spins_examined += 1
                # flip the spin and mark the spin
                self.spin_marks[spin] = self.spin_marker_value
                self.state[spin] = -self.state[spin]

                # visit all links connecting the spin to its nearest neighbors 
                # the bond is activated with probability "probability" if the spin
                # is not equal to its neighbor
                neighbor_spin = (spin + 1) & self.Nm1
                if self.state[spin] != self.state[neighbor_spin] and self.random_nums.next() < probability:
                    self.wolff_stack.push(neighbor_spin)

                neighbor_spin = (spin - 1) & self.Nm1
                if self.state[spin] != self.state[neighbor_spin] and self.random_nums.next() < probability:
                    self.wolff_stack.push(neighbor_spin)

                neighbor_spin = (spin + self.x_size) & self.Nm1
                if self.state[spin] != self.state[neighbor_spin] and self.random_nums.next() < probability:
                    self.wolff_stack.push(neighbor_spin)

                neighbor_spin = (spin - self.x_size) & self.Nm1
                if self.state[spin] != self.state[neighbor_spin] and self.random_nums.next() < probability:
                    self.wolff_stack.push(neighbor_spin)
        
        return spins_examined
    
    cpdef double wolff_sweep_total(self):
        cdef long spins_examined = 0
        cdef long steps = 0
        while spins_examined < self.N:
            spins_examined += self.wolff_sweep()
            steps += 1

        return (<double> spins_examined) /(<double> steps)

    cpdef long magnetization(self):
        """
        calculates the sum of all the spins
        """
        cdef int M = 0
        cdef int i
        for i in range(self.N):
            M += self.state[i]
        
        return M


