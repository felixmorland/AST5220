import matplotlib.pyplot as plt
import numpy as np

x, Theta0, Theta1, Theta2, Phi, delta_CDM = np.loadtxt('perturbation_data/perturbations_k0.01.txt', unpack=True)

plt.semilogy(x, delta_CDM)
plt.show()