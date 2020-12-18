#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np

depth, particles, pions, electrons, photons, muons, neutrinos, energy = np.loadtxt('result.txt', unpack=True)

plt.plot(depth, pions, label='pions')
plt.plot(depth, muons, label='muons')
plt.plot(depth, electrons, label='electrons')
plt.plot(depth, photons, label='photons')
plt.plot(depth, neutrinos, label='neutrinos')
plt.plot(depth, energy / 1e9, label='energy (GeV)')

plt.yscale('log')
plt.legend()
plt.show()