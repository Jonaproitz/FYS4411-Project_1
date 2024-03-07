import numpy as np
import matplotlib.pyplot as plt

filename = "ExpectationValues.txt"
alpha, energy, variance = np.loadtxt(filename, unpack=True)

plt.plot(alpha, variance)
plt.savefig("doc/figures/Varience_B.png")
plt.close()