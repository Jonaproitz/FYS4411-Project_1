import numpy as np
import matplotlib.pyplot as plt

filename = "Results/ExpectationValues_problem_b.txt"
alpha, energy, variance = np.loadtxt(filename, unpack=True)

plt.figure()
plt.title("Variance as a function of parameter $\\alpha$")
plt.plot(alpha, variance)
plt.xlabel("$\\alpha$")
plt.ylabel("Variance")
plt.grid()
plt.savefig("doc/figures/Varience_B.png")
plt.close()