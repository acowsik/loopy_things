from mcmc.mcmc import Ising
import numpy as np
import tqdm
import matplotlib.pyplot as plt

temperature = 2

spins = Ising(10, 10, temperature)

Ms = []

for _ in tqdm.tqdm(range(1024)):
    Ms.append(spins.magnetization() / spins.N)
    spins.sweep()

plt.plot(Ms)
plt.figure()
plt.plot(np.correlate(Ms, Ms, mode='same'))
plt.show()