import numpy as np
import matplotlib.pyplot as plt

numberOfParticles = 1
filename = f"Results/ExpectationValues_problem_b_{numberOfParticles}.txt"
alpha, energy, variance = np.loadtxt(filename, unpack=True)

plt.rcParams.update({'font.size': 12})

plt.figure()
plt.title(f"Variance for {numberOfParticles} particles in 3 dimensions")
plt.plot(alpha, variance)
plt.xlabel("$\\alpha$")
plt.ylabel("Variance")
plt.grid()
plt.savefig(f"doc/figures/Varience_B_{numberOfParticles}.png")
plt.close()

plt.figure()
plt.title(f"Energy for {numberOfParticles} particles in 3 dimensions")
plt.plot(alpha, energy)
plt.xlabel("$\\alpha$")
plt.ylabel("Energy")
plt.grid()
plt.savefig(f"doc/figures/Energy_B_{numberOfParticles}.png")
plt.close()