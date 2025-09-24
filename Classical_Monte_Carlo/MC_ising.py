import numpy as np

# Parameters 
class Ising_MC_2D:
    def __init__(self, Nx, Ny, Temp, nwarm=800, nmeas=500):
        # Initialize constants 
        self.Nx = Nx            # Width of lattice
        self.Ny = Ny            # Height of lattice
        self.N = Nx*Ny          # Number of sites
        self.T = Temp           # Temperature
        self.nwarm = nwarm      # Number of warm-up to reach thermal equilibrium
        self.nmeas = nmeas      # Number of Monte Carlo measurments 

        # initialize lookup table
        self.prob_table= {4: np.exp(-8/self.T),
                          2: np.exp(-4/self.T)} 

        # Initiate spin configuration to be updated during the simulation
        self.spins = np.ones((Nx,Ny))
        #self.spins = np.random.choice([-1, 1], size=(Nx,Ny))

        # Initialize measurments 
        self.magnetization = 0
        self.magnetization_squared = 0 
        self.energy = 0
        self.energy_squared = 0

        # Initialize energy and update it during the simulation to reduce repeated calculations
        self.energy_help = 0
        for i in range(Nx):
            for j in range(Ny):
                self.energy_help -= 0.5*self.spins[i,j]*(self.spins[(i+1)%Nx,j] +self.spins[(i-1)%Nx,j] + self.spins[i, (j-1)%Ny] +self.spins[i, (j+1)%Ny])

    #def calculate_energy(self):
    #    energy = 0
    #    for i in range(self.Nx):
    #        for j in range(self.Ny):
    #            self.energy -= 0.5*self.spins[i,j]*(self.spins[(i+1)%self.Nx,j] + self.spins[(i-1)%self.Nx,j]+ self.spins[i, (j-1)%self.Ny] +self.spins[i, (j+1)%self.Ny])               
    #    return energy

    def run_MC_simulation(self):
        for sweep in range(self.nwarm + self.nmeas):
            for i in range(self.Nx):
                for j in range(self.Ny):
                    # Instead picking random sites we do full lattice sweeps 
                    # Calculate energy difference 
                    h = self.spins[(i+1)%self.Nx, j] + self.spins[(i-1)%self.Nx, j] + self.spins[i, (j+1)%self.Ny] + self.spins[i, (j-1)%self.Ny] 
                    if h*self.spins[i,j] <= 0:
                        self.energy_help += 2*h*self.spins[i,j]
                        self.spins[i,j] *= -1
                    else:
                        p = np.random.rand()
                        if p <= self.prob_table.get(h*self.spins[i,j]):
                            self.energy_help += 2*h*self.spins[i,j]
                            self.spins[i,j] *= -1
            # Perform measurments 
            if sweep > self.nwarm:
                # Magentic measuremnts 
                mag = np.sum(self.spins)/self.N
                self.magnetization += mag
                self.magnetization_squared += mag**2 

                # Energy Measurments 
                self.energy += self.energy_help
                self.energy_squared += (self.energy_help)**2
        



