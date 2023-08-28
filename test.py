from mcmc.mcmc import Ising, Stack
import numpy as np
import tqdm
import matplotlib.pyplot as plt

s = Stack(20)
s.push(1)
s.push(4)
s.push(100)

print(s.pop())
print(s.pop())
print(s.pop())

temperatures = np.linspace(2.2, 2.3, 10)
magnetizations = []
mag_squareds = []
for temperature in temperatures:
    spins = Ising(10, 10, temperature)
    for _ in range(20):
        spins.wolff_sweep_total()

    Ms = []

    for i in tqdm.tqdm(range(1024)):
        if i % 32 == 0:
            Ms.append(spins.magnetization() / spins.N)
        spins.wolff_sweep_total()
    
    magnetizations.append(np.mean(Ms))
    mag_squareds.append(np.std(Ms))

plt.figure()
plt.plot(temperatures, magnetizations)
plt.ylabel("M")
plt.figure()
plt.plot(temperatures, mag_squareds)
plt.ylabel("$M^2$")
plt.show()