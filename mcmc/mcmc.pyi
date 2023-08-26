# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
from typing import List
import numpy as np

class Ising:
    state: np.ndarray
    x_size: int
    y_size: int
    N: int
    beta: float

    def __init__(self, log_x_size: int, log_y_size: int, temperature: float) -> None: ...
        # Initialize the state with x_size 2**log_x_size and y_size 2**log_y_size and at the given temperature 
    def step(self) -> None: ...
        # does one potential spin flip 
    def sweep(self) -> None: ...
        # does N potential spin flips
    def magnetization(self) -> int: ...
        # calculates the sum of all the spins


