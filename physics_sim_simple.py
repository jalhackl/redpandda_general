import itertools
import numpy as np

timestep = 0.01
hconst = 5
hconst2 = 2

dyn_friction_on = 0

partdim = 3



simplot = 0

#kt = 1.38e-23
kt = 1

#arbitrarily chosen constants for the potential
hconst = 20
hconst2 = 0.5

#turn harmonic potential on (1) or off (0)
harmpot_on=1

#turn fene potential on (1) or off (0)
fenepot_on=0

#turn quartic potential on (1) or off (0)
quarticpot_on=0

#turn thermostat on (1) or off (0)
thermostat_on = 0

# space-dependent friction
dyn_friction_on = 0

#constants for friction term 
frica = 0.5
fricb = 5

feneK = hconst


class Particle():
    newid = itertools.count().__next__
    def __init__(self, position, velocity, mass, friction , force , temperature):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.friction = friction
        self.force = force
        self.temperature = temperature
        
        self.positions = []
        self.velocities = []
        self.forces = []        
        self.kinetic1 = []
        self.potentiale = []
        self.frictions = []
        
        
        
        # additional friction for specific particles when using dynamic friction (can be assigned for each particle)
        self.frictionaadd = 2
        self.frictionbadd = 2
        self.frictionamult = 1
        self.frictionbmult= 1
        
        #friction at initialization
        self.startfriction = friction
        
        self.id = Particle.newid()
        
        #variables for 2nd order integrator
        self.c=0
        self.theta=0
        self.xi=0
        
        #counting of steps performed
        self.counter = 0


    
#standard velocity verlet scheme: new positions
def new_pos(Particle):
    
    if (dyn_friction_on == 1):
        #dynamic_frictions(Particle)
        Particle.counter = Particle.counter + 1
        if ((Particle.counter % 1000000) == 0):
            print("counter:")
            print(Particle.counter)
            print("Friction:")
            print(Particle.friction)
            print("position:")
            print(Particle.position)
            print("velocity:")
            print(Particle.velocity)
       
    Particle.positions.append(Particle.position) 
    Particle.velocities.append(Particle.velocity) 

    Particle.position = Particle.position + Particle.velocity*timestep + (Particle.force )*(timestep**2)*0.5 / Particle.mass


#standard velocity verlet scheme: new velocities
def new_vel(Particle):

    actforces = new_forces_thermostat(Particle)  

    Particle.velocity = Particle.velocity + timestep * (actforces + Particle.force) / (2*Particle.mass)
    
    Particle.force = actforces
    

#updating of the acting forces
def new_forces_thermostat(Particle):
    force = 0
    force = force #- Particle.friction * Particle.velocity * Particle.mass
    if (harmpot_on == 1):
        hp = harm_potential(Particle)
        force = force - hp
        
    if (fenepot_on == 1):
        fp = fene_potential(Particle)
        force = force - fp
       # print("fenep")
       # print(force)
    if (quarticpot_on == 1):
        qp = quartic_potential(Particle)
        force = force - qp
    
    force = force + np.random.randn(partdim)*sigma_thermostat(Particle)
    
    Particle.forces.append(Particle.force)
   
    return force


#Langevin thermostat
#see https://www.vasp.at/wiki/index.php/Langevin_thermostat
#this function computes the variance of the random number used for the random force in the Langevin equation
#via the 'temperature' attribute the temperature can be adjusted for each particle
def sigma_thermostat(Particle):
    if (thermostat_on == 1):
        sigma = np.sqrt(2*kt*Particle.temperature*Particle.friction*Particle.mass/timestep)
    else:
        sigma = 0
    return sigma


#the potential between two particles
def harm_potential(Particle):
    pot = 0
    pote = 0
    for ParticleP in particles:
        
        if (id(Particle) == id(ParticleP)):
            
            continue
        else:          
            dst = ParticleP.position-Particle.position
           
            vk = dst / np.linalg.norm(dst)
            
            pott = hconst * 2 * np.linalg.norm(dst - hconst2)
            pot = pot - pott * vk
     
            pot = pot / 2

    return pot



#the potential between two particles
def harm_potential(Particle):
    pot = 0
    pote = 0
    for ParticleP in particles:
        
        if (id(Particle) == id(ParticleP)):
            
            continue
        else:          
            #dst = minimagedist(Particle.position, ParticleP.position, "p")
            dst = ParticleP.position-Particle.position
           
            vk = dst / np.linalg.norm(dst)
            
            #pott = hconst * 2 * np.linalg.norm(dst - hconst2)
            #pott = hconst * np.linalg.norm(dst - hconst2)
            
            pott = hconst * (np.linalg.norm(dst) - hconst2)

            pot = pot - pott * vk
     
            pot = pot # / 2

    return pot


#fene potential
def fene_potential(Particle):
    pot = 0
    
    K = feneK
    R0 = 1.5
    epsilon = 1
    sigma_fene = 1
    cutoff = 2 **(1/6)
    #fene potential (https://lammps.sandia.gov/doc/bond_fene.html#examples)
    # E = -0.5 * K * R0^2 * ln(1-(r/R0)^2)) + LJ + epsilon
    # -> force:
    # F = K * r / (1 - (r / R0)^2) + (24 * epsilon / r) * (2 (sigma_fene / r)^12 - (sigma_fene / r)^6)
    
    for ParticleP in particles:
        if (id(Particle) == id(ParticleP)):
            
            continue
        else: 
            dst = ParticleP.position-Particle.position
            r = np.linalg.norm(dst)
            
            vk = dst / np.linalg.norm(dst)
            if (r < cutoff):
                pott = K * (1 / (1 - (r / R0)**2)) - (24 * epsilon / r) * (2 * (sigma_fene / r)**12 - (sigma_fene / r)**6) 
            else:
                pott = K * (1 / (1 - (r / R0)**2)) 
           
            pot = pot - pott * vk 
   
    return pot



#quartic potential
def quartic_potential(Particle):
    pot = 0
    for ParticleP in particles:
        if (id(Particle) == id(ParticleP)):
            
            continue
        else: 
            dst = ParticleP.position-Particle.position
            r = np.linalg.norm(dst)
            vk = dst / np.linalg.norm(dst)
            
            pott = (hconst * 4 * r**3 ) 
            
            pot = pot - pott * vk 
    return pot



def run_simulation(particles):
    
    # loop for 2 (or more) particles  
    for step in range(steps):
        
        
        for Particle in particles:        
            # new positions
            new_pos(Particle)
        for Particle in particles:  
            # new velocities
            new_vel(Particle)
        


#simulation parameters: steps and stepsize
steps = 1000
dt = 0.01
timestep = dt




def init_particles(temp2 = 5):
    global particle5
    global particle6
    global partdim
    global numpart
    global particles
    #initialization of particles:
    #arguments: 
    # position
    # velocity
    # mass
    # friction
    # intial force acting on particle (set to 0)
    # temperature (which is enforced by the thermostat, if switched on)

    # set friction coefficient for all particles    
    setfriction = 1

    #1d
    #particle5 = Particle(np.array([1]),np.array([0]),1,np.array([setfriction]),np.array([0]),1)
    #particle6 = Particle(np.array([2]),np.array([0]),1,np.array([1]),np.array([0]),temp2)


    #3d
    particle5 = Particle(np.array([1,2,-1]),np.array([0,0,0]),1,np.array([1,1,1]),np.array([0,0,0]),1)
    particle6 = Particle(np.array([1.75,2.5,-1.5]),np.array([0,0,0]),1,np.array([setfriction,setfriction,setfriction]),np.array([0,0,0]),temp2)

    partdim = len(particle5.position)
    
    #add all particles to a set
    particles = set()
    particles.add(particle5)
    particles.add(particle6)
    numpart = len(particles)
    
    return particles









temp2 = 5
particles = init_particles(temp2)
    

run_simulation(particles)

