#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np

depth, particles, nuclei, pions, electrons, photons, muons, neutrinos = np.loadtxt('result.txt', unpack=True)

plt.plot(depth, nuclei, label='nuclei')
plt.plot(depth, pions, label='pions')
plt.plot(depth, muons, label='muons')
plt.plot(depth, electrons, label='electrons')
plt.plot(depth, photons, label='photons')
plt.plot(depth, neutrinos, label='neutrinos')

plt.yscale('log')
plt.legend()
plt.show()