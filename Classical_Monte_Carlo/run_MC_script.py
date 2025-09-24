from MC_ising import * 
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


#ising = Ising_MC_2D(20,20,1, 1000, 1000)
#ising.run_MC_simulation()
#print('Done')

#print(rf"Energy: {ising.energy/ising.nmeas}")
#print(rf"Magnetization {ising.magnetization/ising.nmeas}")

Nx = 20
Ny = 20
N = Nx*Ny
T_vals = np.linspace(1,4,100)
magnitization = np.zeros_like(T_vals)
mag_sus = np.zeros_like(T_vals)
energies = np.zeros_like(T_vals)
heat_capacity = np.zeros_like(T_vals)

magnitization_err = np.zeros_like(T_vals)
mag_sus_err = np.zeros_like(T_vals)
energies_err = np.zeros_like(T_vals)
heat_capacity_err = np.zeros_like(T_vals)


for i, temp in enumerate(tqdm(T_vals)):
    mag_ = []
    energy_ = []
    chi_ = []
    cv_ = []
    for _ in range(15): #independent markov chains
        MC_simulation = Ising_MC_2D(Nx, Ny, temp, 20, 700)
        MC_simulation.run_MC_simulation()
        mag_.append(abs(MC_simulation.magnetization/MC_simulation.nmeas))
        chi_.append(N*(MC_simulation.magnetization_squared/MC_simulation.nmeas - (MC_simulation.magnetization/MC_simulation.nmeas)**2)/temp)
        energy_.append(MC_simulation.energy/MC_simulation.nmeas)
        cv_.append((MC_simulation.energy_squared/MC_simulation.nmeas - (MC_simulation.energy/MC_simulation.nmeas)**2)/(temp*temp))

    # Magnetic properties 
    magnitization[i] += np.mean(mag_)
    magnitization_err[i] += np.sqrt(np.var(mag_))

    mag_sus[i] += np.mean(chi_)
    mag_sus_err[i] += np.sqrt(np.var(chi_))
    # Energy properties 
    energies[i] += np.mean(energy_)
    heat_capacity[i] += np.mean(cv_)


# Make plots
fig, axs = plt.subplots(2,2, figsize=(10,8))

# Top left
onsager_mag = (1-(np.sinh(2/T_vals))**(-4))**(1/8)
axs[0,0].plot(T_vals, abs(magnitization), '.')
axs[0,0].plot(T_vals, onsager_mag, color='red', alpha=0.5)
axs[0,0].set_ylabel(r"$M(T)$")
axs[0,0].set_xlim(1,4)

# Top right
axs[0,1].plot(T_vals, energies/N, '.')
axs[0,1].set_ylabel(r"$E/N$")
axs[0,1].set_xlim(1,4)

# Bottom left
axs[1,0].plot(T_vals, mag_sus, '.')
axs[1,0].set_ylabel(r"$\chi(T)$")
axs[1,0].set_xlabel(r"$T$")
axs[1,0].set_xlim(1,4)


axs[1,1].plot(T_vals, heat_capacity/N, '.')
axs[1,1].set_ylabel(r"$C_V(T)$")
axs[1,1].set_xlabel(r"$T$")
axs[1,1].set_xlim(1,4)

plt.show()
