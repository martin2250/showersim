#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np

depth, particles, nuclei, pions, electrons, photons, muons, neutrinos, remainin_energy, ionization = np.loadtxt('result.txt', unpack=True)

plt.plot(depth, nuclei, label='nuclei')
plt.plot(depth, pions, label='pions')
plt.plot(depth, muons, label='muons')
plt.plot(depth, electrons, label='electrons')
plt.plot(depth, photons, label='photons')
plt.plot(depth, neutrinos, label='neutrinos')
plt.plot(depth, remainin_energy / 1e9, label='remainin_energy (GeV)')
plt.plot(depth, ionization / 1e9, label='ionization (GeV/m)')

plt.yscale('log')
plt.legend()
plt.show()