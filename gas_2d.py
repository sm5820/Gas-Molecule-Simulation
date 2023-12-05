#Object-oriented implementation of a hard-disks molecular dynamics simulations
from itertools import combinations

class Vector():
    '''2D vectors'''

    def __init__(self, i1,i2):
        '''Initialise vectors with x and y coordinates'''
        self.x = i1
        self.y = i2

    def __add__(self,other):
        '''Use + sign to implement vector addition'''
        return (Vector(self.x+other.x,self.y+other.y))

    def __sub__(self,other):
        '''Use - sign to implement vector "subtraction"'''
        return (Vector(self.x-other.x,self.y-other.y))

    def __mul__(self,number):
        '''Use * sign to multiply a vector by a scaler on the left'''
        return Vector(self.x*number,self.y*number)

    def __rmul__(self,number):
        '''Use * sign to multiply a vector by a scaler on the right'''
        return Vector(self.x*number,self.y*number)

    def __truediv__(self,number):
        '''Use / to multiply a vector by the inverse of a number'''
        return Vector(self.x/number,self.y/number)

    def __repr__(self):
        '''Represent a vector by a string of 2 coordinates separated by a space'''
        return '{x} {y}'.format(x=self.x, y=self.y)

    def copy(self):
        '''Create a new object which is a copy of the current.'''
        return Vector(self.x,self.y)

    def dot(self,other):
        '''Calculate the dot product between two 2D vectors'''
        return self.x*other.x + self.y*other.y

    def norm(self):
        '''Calculate the norm of the 2D vector'''
        return (self.x**2+self.y**2)**0.5

class Particle():
    def __init__(self, position, momentum, radius, mass):
        'Initialises new particle by assigning position, momentum, radius and mass' 
        self.position = position 
        self.momentum = momentum 
        self.radius = radius 
        self.mass = mass 
        
    def velocity(self):
        return (self.momentum/self.mass)
    
    def kinetic_energy(self):
        return((0.5*self.momentum.dot(self.momentum))/self.mass)
    
    def copy(self):
        return Particle(self.position, self.momentum, self.radius, self.mass)
    
    def overlap(self, other_particle):
        distance = (self.position - other_particle.position).norm()
        radii_sum = self.radius + other_particle.radius
        return distance < radii_sum

class Simulation():
    def __init__(self, particles, box_length, dt):
        'Initialises ... ' 
        self.particles = particles
        self.box_length = box_length
        self.dt = dt
        self.trajectory = []
    
    def apply_particle_collisions(self, particle1, particle2):
        momentum1 = particle1.momentum.copy()
        momentum2 = particle2.momentum.copy()
        mass1 = particle1.mass
        mass2 = particle2.mass
        momentum = momentum2-momentum1
        pos_vec = particle2.position - particle1.position 
        if momentum.dot(pos_vec)<0:
            if particle1.overlap(particle2):
                norm_pos_vec = pos_vec/pos_vec.norm()
            
            
            
                collision_comp1 = norm_pos_vec * momentum1.dot(norm_pos_vec)
                collision_comp2 = norm_pos_vec * momentum2.dot(norm_pos_vec)
            
                ortho1 = momentum1 - collision_comp1
                ortho2 = momentum2 - collision_comp2
            
                collision_comp1_copy = collision_comp1
                collision_comp2_copy = collision_comp2
            
                collision_comp1 = ((mass2-mass1)*collision_comp1 + 2*mass1*collision_comp2_copy)/(mass1+mass2)
                collision_comp2 = ((mass2-mass1)*collision_comp2 + 2*mass2*collision_comp1_copy)/(mass1+mass2)
            
                particle1.momentum = ortho1 + collision_comp1
                particle2.momentum = ortho2 + collision_comp2
            
    
    
    def apply_box_collisions(self, particles):
        for particle in self.particles:
            if (particle.position.x - particle.radius) <= 0 and particle.momentum.x <0:
                particle.position.x = particle.radius
                particle.momentum.x = abs(particle.momentum.x)
            if (particle.position.y- particle.radius) <= 0 and particle.momentum.y <0:
                particle.momentum.y = abs(particle.momentum.y)
                particle.position.y = particle.radius
            if (particle.position.y + particle.radius) >= self.box_length and particle.momentum.y >0:
                particle.position.y = self.box_length - particle.radius
                particle.momentum.y= -abs(particle.momentum.y)
            if (particle.position.x + particle.radius) >= self.box_length and particle.momentum.x > 0:
                particle.position.x = self.box_length - particle.radius
                particle.momentum.x=  -abs(particle.momentum.x)
    
    def step(self):
        for particle in self.particles:
            particle.position = particle.position + particle.velocity()*self.dt
            self.apply_box_collisions(self.particles)
        for particle1, particle2 in combinations(self.particles, 2): 
            self.apply_particle_collisions(particle1,particle2)
        self.record_state()
                
    def record_state(self):
        state = []
        for particle in self.particles:
            state.append(particle.copy())
        self.trajectory.append(state)
