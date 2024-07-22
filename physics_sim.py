import numpy as np

simplot = 0
simplot3d = 0

#kt = 1.38e-23
kt = 1

#arbitrarily chosen constants for the potential
hconst = 30
hconst2 = 1.5

#turn harmonic potential on (1) or off (0)
harmpot_on=1

#turn fene potential on (1) or off (0)
fenepot_on=0

#turn quartic potential on (1) or off (0)
quarticpot_on=0

#turn thermostat on (1) or off (0)
thermostat_on = 1
#include friction thermostat, but not Langevin thermostat 
onlyfriction = 0

# space-dependent friction
dyn_friction_on = 0

#constants for friction term 
frica = 0.5
fricb = 5

feneK = hconst

#hydrodynamics
sph = 0
mpcd = 0

#box = ((-1e-9mp, 1e-9), (-1e-9, 1e-9))

#WALDER 15
boxhl = 5 
#boxhl = 5


boxa = -boxhl
boxe = boxhl
    
box = ((-boxhl, boxhl), (-boxhl, boxhl), (-boxhl, boxhl))

boxl = np.abs(boxhl) + np.abs(boxhl)
#boxcheck-unterarten

boxch = 1
boxchsimple = 0
boxchperiodiccentered = 1
boxchflow = 0


#
nearconv = 1

thermocollision =0
langevinthermocollision = 0

#naechstes egal fuer rotation matrix
collisionuseparticletemp = 0

usealwayssolvt = 1
useavgcelltemp = 0

rotcellscaletemp = 1
alphaangle = 130


leesedwards = 0

shiftparticlestoorigin = 0

followboxes = 1
if (followboxes == 1):
    boxesinfo = []
boxcounter = 0


fluidinfo = []

#numberboxes = 3375
mpcd_boxl = 1.0
numberboxes = (boxl / mpcd_boxl ) ** 3

#double length
#numberboxes = 27000


#numberboxes = 4096
#numberboxes = 64


diffcoeffs = []

solvent_energies = []

#fuer polymer
solvent_energies_poly = []

shearnumbers_array = []

invariance_shift = True

boxes1d = numberboxes ** (1./3.)
mpcdboxl = boxl / boxes1d

randomv = np.random.uniform(-mpcdboxl/2, mpcdboxl/2, 3)

polymernr = 0

check_momentum = 1
check_momentum_cell = 1

total_momentum_pre = []
total_momentum = []

impl_obstacles = False

SABPO= False

EABPO = True

LAN = False

alan = False

changeek = True
EABPO_change_fluid = True
EABPO_mono_backflow = True
EABPO_mono_active_backflow = True
#EABPO_fluid_backflow = False



class Particle():
    newid = itertools.count().__next__
    def __init__(self, position, velocity, mass, friction , force , temperature, solvent = 0, mpcd_off=False, const_vel=False, active_velocity0 = 1, active_velocity_vec=np.array([1,1,1]), sab=False, eab=False, lan=False, special=False, ala=False, thermostatted=False):
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
        
        self.sab = sab
        self.eab = eab
        
        self.lan = lan
        
        self.alan = ala
        
        self.mpcd_off = mpcd_off
        self.const_vel = const_vel
        
        self.lantemperature = False
        
        self.thermostatted = thermostatted
        
        #SABPO
        if (SABPO == True and self.sab == True) or alan==True:
            #v0
            self.active_velocity0 = active_velocity0
            #ek
            self.active_velocity_vec = active_velocity_vec / np.linalg.norm(active_velocity_vec)
            #active_velocity is given by v0*ek
            self.active_velocity = self.active_velocity0 * self.active_velocity_vec
            
            self.active_velocities = []
            self.active_velocities_vec = []
            self.total_velocities = []
            
            self.total_velocity = self.velocity + self.active_velocity
            
        if (EABPO == True and self.eab == True) or alan==True:
            #v0
            self.active_velocity0 = active_velocity0
            #ek
            self.active_velocity_vec = active_velocity_vec / np.linalg.norm(active_velocity_vec)
            #active_velocity is given by v0*ek
            #self.active_velocity = self.active_velocity0 * self.active_velocity_vec
            
            self.active_force = self.friction[0] * active_velocity0 * active_velocity_vec
            
            #self.active_velocities = []
            self.active_velocities_vec = []
            #self.total_velocities = []
            
            self.active_forces = []
            #self.total_velocity = velocity + self.active_velocity
            
        if (EABPO == True or SABPO == True) and self.sab == False and self.eab == False:
            self.active_velocity0 = 0
            self.active_velocity = np.array([0.0,0.0,0.0])
            self.active_force = np.array([0.0,0.0,0.0])
            self.total_velocity = self.velocity + self.active_velocity
            
            self.active_velocities = []
            self.total_velocities = []
            
            self.active_velocity_vec = active_velocity_vec
            #self.active_velocity_vec = np.zeros(3)
            
            self.active_velocities_vec = []
            
            self.active_forces = []
            
            
            #Particle.active_velocities_vec.append(Particle.active_velocity_vec)
        if (EABPO == True or SABPO == True) and (self.sab == True or self.eab == True) and self.lan == True:
            self.temperature = self.active_velocity0
            
            self.lantemperature = self.active_velocity0
            
            self.active_force = np.array([0.0,0.0,0.0])
            self.active_forces = []
            
        
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
        
        #hydro
        #if solvent == 1:
        #    self.solvent = solvent
        #else:
        #    self.solvent = 1
        self.solvent = solvent
            
        
        if solvent == 1:
            self.polynr = 0
        else:
            self.polynr = polymernr
        
        
        
        #counting of steps performed
        self.counter = 0
        
        #for computation of msd (unbound positions)
        self.boxcross = [0,0,0]
        
        self.startposition = position
        
        self.special = special
        
   


   #standard velocity verlet scheme: new positions
def new_pos(Particle):
    
    oldpos = Particle.position
        
    if (dyn_friction_on == 1):
        dynamic_frictions(Particle)
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
    
    if (Particle.solvent == 0):
        Particle.positions.append(Particle.position) 
        Particle.velocities.append(Particle.velocity) 
        
        

        
        if (Particle.solvent == 0 and (SABPO == True or Particle.alan == True)):
            Particle.active_velocities.append(Particle.active_velocity)
            Particle.total_velocity = Particle.velocity + Particle.active_velocity
            Particle.total_velocities.append(Particle.total_velocity)
        

        
       

    if (Particle.solvent == 0):
        Particle.position = Particle.position + Particle.velocity*timestep + (Particle.force )*(timestep**2)*0.5 / Particle.mass
    
        #alan
        if Particle.alan == True:
            Particle.position = Particle.position + Particle.active_velocity*timestep + (Particle.active_force )*(timestep**2)*0.5 / Particle.mass

            
        
        #SABPO
        if SABPO == True and Particle.sab == True and Particle.lan == False:
            Particle.position = Particle.position + timestep*Particle.active_velocity
          
        if SABPO == True and Particle.sab == True and Particle.lan == True:
            Particle.position = Particle.position + (Particle.active_velocity )*(timestep**2)*0.5 / Particle.mass # timestep*Particle.active_velocity
        
        if EABPO == True and Particle.eab == True and Particle.mpcd_off == False:
           # if changeek == True:
           #     Particle.active_velocity_vec = integrate_evel_vector(Particle, Particle.active_velocity_vec)
            
           # Particle.active_force = Particle.friction[0] * active_velocity0 * active_velocity_vec
            
            Particle.position = Particle.position + (Particle.active_forces[-1] )*(timestep**2)*0.5 / Particle.mass
            
           # global massparts
           # global total_active_force 
            
            if (EABPO_mono_active_backflow == True):
                backflow_force = - (Particle.mass/massparts) * old_total_active_force

                Particle.position = Particle.position + (backflow_force)*(timestep**2)*0.5 / Particle.mass
            
        if EABPO == True and Particle.eab == False and EABPO_mono_backflow == True and Particle.mpcd_off == False:
            backflow_force = - (Particle.mass/massparts) * old_total_active_force

            #Particle.position = Particle.position - ((timestep ** 2)/(2*massparts)) * total_active_force
            Particle.position = Particle.position + (backflow_force)*(timestep**2)*0.5 / Particle.mass
    
    else:
        Particle.position = Particle.position + Particle.velocity*mpcdtimestep #+ (Particle.force )*(mpcdtimestep**2)*0.5 / Particle.mass

        if (EABPO == True and EABPO_change_fluid == True and Particle.mpcd_off == False):
           # global massparts
           # global total_active_force 
            
            #Particle.position = Particle.position - ((mpcdtimestep ** 2)/(2*massparts)) * mpcd_total_active_force
            backflow_force = - (Particle.mass/massparts) * old_total_active_force

            #Particle.position = Particle.position - ((mpcdtimestep ** 2)/(2*massparts)) * mpcd_total_active_force
            
            Particle.position = Particle.position + (backflow_force)*(mpcdtimestep**2)*0.5 / Particle.mass
            
        
 #new
    if (Particle.solvent == 0 and impl_obstacles == True):

        impl_obstacles_calc(Particle, oldpos)   



#EABPO
def compute_active_forces(Particle):
    if changeek == True and Particle.lan == False:
                         
        Particle.active_velocities_vec.append(Particle.active_velocity_vec)

        
        Particle.active_velocity_vec = integrate_evel_vector(Particle, Particle.active_velocity_vec)
        
    if Particle.lan == True:
            new_temp = new_forces_thermostat_lan(Particle) 
            
            
    
    Particle.active_forces.append(Particle.active_force)
              
    if Particle.lan == True:
        Particle.active_force = new_temp
    else:
        Particle.active_force = Particle.friction[0] * Particle.active_velocity0 * Particle.active_velocity_vec


#EABPO
def get_backflow_force(particlespoly):
    global total_active_force
    
    total_active_force = np.array([0.0,0.0,0.0])
    for Particle in particlespoly:
        if Particle.mpcd_off == False:
            total_active_force = total_active_force + Particle.active_force
        
    return total_active_force



#standard velocity verlet scheme: new velocities
def new_vel(Particle):
    if (Particle.solvent == 0):
        actforces = new_forces_thermostat(Particle) 
        if (sph ==1):
            actforces = actforces + getAcc(Particle)
    else:
        actforces = 0
        Particle.force = actforces
        if (EABPO == True and EABPO_change_fluid == True):
            global massparts
            global total_active_force 
            
            #Particle.velocity = Particle.velocity - (mpcdtimestep / massparts) * mpcd_total_active_force
            
            Particle.velocity = Particle.velocity - mpcdtimestep * ( (Particle.mass / massparts) * ( total_active_force + old_total_active_force ) ) / (2*Particle.mass)

            return
            
        
        else:
            return
        

        if (sph ==1):
            actforces = actforces + getAcc(Particle)

        
    Particle.velocity = Particle.velocity + timestep * (actforces + Particle.force) / (2*Particle.mass)
    
    #alan
    if Particle.alan == True:
        
        new_temp = new_forces_thermostat_lan(Particle) 
        Particle.active_force = new_temp
        Particle.active_forces.append(new_temp)
        Particle.active_velocity + timestep * (Particle.active_force + Particle.active_forces[-1]) / (2*Particle.mass)
    
    #EABPO: backflow force on non-active polymer particles
    if (EABPO == True and Particle.eab == True and Particle.mpcd_off == False):
        Particle.velocity = Particle.velocity + timestep * (Particle.active_force + Particle.active_forces[-1]) / (2*Particle.mass)

    
    if (EABPO == True):
        if ((Particle.eab == True and EABPO_mono_active_backflow == True and Particle.mpcd_off == False) or (Particle.eab == False and EABPO_mono_backflow == True and Particle.mpcd_off == False)):
            Particle.velocity = Particle.velocity - timestep * ( (Particle.mass / massparts) * ( total_active_force + old_total_active_force ) ) / (2*Particle.mass)


    #SABPO
    if (Particle.solvent == 0 and SABPO == True and Particle.sab == True):
        
        
        Particle.active_velocities_vec.append(Particle.active_velocity_vec)
        
        
        if Particle.lan == False:
            new_evel = integrate_evel_vector(Particle, Particle.active_velocity_vec)
            new_active_velocity = Particle.active_velocity0 * new_evel

            Particle.active_velocity_vec = new_evel

            #Particle.active_velocity_vec.append(new_evel)
        else:
            new_temp = new_forces_thermostat_lan(Particle) 
            
            
            
            new_active_velocity = Particle.active_velocity + timestep * (new_temp + Particle.active_force)  / (2*Particle.mass)
            
            Particle.active_force = new_temp
            
            
            Particle.active_forces.append(new_temp)

        
        Particle.active_velocity = new_active_velocity
        
        #2212 ausk
        #Particle.total_velocity = Particle.velocity + new_active_velocity
        #Particle.total_velocities.append(Particle.total_velocity)
    

    if const_vel_try == True:
        Particle.velocity = const_vel
        
    Particle.force = actforces
    


