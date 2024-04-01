import numpy as np
import matplotlib.pyplot as plt

numberOfParticles = 1
filename = f"Results/ExpectationValues_problem_b_{numberOfParticles}.txt"
alpha, energy, variance = np.loadtxt(filename, unpack=True)

plt.figure()
plt.title("Variance as a function of parameter $\\alpha$")
plt.plot(alpha, variance)
plt.xlabel("$\\alpha$")
plt.ylabel("Variance")
plt.grid()
plt.savefig(f"doc/figures/Varience_B_{numberOfParticles}.png")
plt.close()