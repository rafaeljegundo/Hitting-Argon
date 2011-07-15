#!/usr/bin/env python
# encoding: utf-8

from random import random
from numpy import pi, sin, cos, sqrt, hypot, arctan
from time import time


# Constantes Universais

c = 299792458 # m/s
K = 8.61E-5 # eV/K
T = 293 # K
argonMass = 931.46E6/(c**2) # eV
electronMass = 511E3/(c**2) # eV
n = 2.18E25 # àtomos /m^3

class Particle:	
	
	"""Ion attributes: position (x,y), velocity (vx,vy) and energy"""
	
	def __init__(self, particle_type, energy, x, y, z, vx, vy, vz):
		
		self.teta = arctan(vy/vx)
		if vz == 0:
			self.phi = 0	
		else:
			self.phi = arctan(hypot(vx,vy)/vz)
		
		if particle_type == "Argon+":		
			self.mass = argonMass
		elif particle_type == "Electron":
			self.mass = electronMass
		
		self.particle_type = particle_type
		
		self.x = x
		self.y = y
		self.z = z
		
		self.energy = energy
		
		self.v = sqrt((2*energy)/self.mass)
		
		q = cos(self.teta)*sin(self.phi)
		u = sin(self.teta)*sin(self.phi)
		s = cos(self.phi)

		self.vx = self.v*q
		self.vy = self.v*u
		self.vz = self.v*s
				
		self.collisioncounter = 0	
				
		return

	def colides(self):
		
		self.collisioncounter += 1
				
		# Generate target atom of the gas
		vxg = random()
		vyg = random()
		vzg = random()
		
		gas_atom = Particle("Argon+", (3/2)*K*T, self.x, self.y, self.z, vxg, vyg, vzg)
		
		#q_gas = sin(gas_atom.phi)*cos(gas_atom.phi)
		#u_gas = sin(gas_atom.phi)*sin(gas_atom.phi)
		#s_gas = cos(gas_atom.phi)
		
		# Conversão para o sistema do centro de massa
		self.vfx = self.vx - gas_atom.vx
		self.vfy = self.vy - gas_atom.vy
		self.vfz = self.vz - gas_atom.vz

		self.vf = sqrt(self.vfx**2 + self.vfy**2 + self.vfz**2)
				
		q = self.vfx/self.vf
		u = self.vfy/self.vf
		s = self.vfz/self.vf
		
		# Qual o m e M?
		gama2 = self.mass/(self.mass+gas_atom.mass)
		gama1 = gas_atom.mass/(self.mass+gas_atom.mass)

		
		#vcmx = gama2*self.vf*q # mas estes dois ultimos termos == vx !!! 
		#vcmy = gama2*self.vf*u
		#vcmz = gama2*self.vf*s
		
		ro = sqrt(1-s**2)
		
		# 
		phi = gas_atom.phi
		teta = gas_atom.teta
		
		# Cosenos directores ' - colisao
#		ql = (sin(self.phi)*cos(self.teta)*s*q-sin(self.phi)*sin(self.teta)*u)/ (ro + cos(self.phi)*q)
#		ul = (sin(self.phi)*cos(self.teta)*s*u-sin(self.phi)*sin(self.teta)*q)/ (ro + cos(self.phi)*u)
#		sl = -sin(self.phi)*sin(self.teta)*ro + cos(self.phi)*s

		# OR
		
		ql = (sin(phi)*cos(teta)*s*q-sin(phi)*sin(teta)*u)/ (ro + cos(phi)*q)
		ul = (sin(phi)*cos(teta)*s*u-sin(phi)*sin(teta)*q)/ (ro + cos(phi)*u)
		sl = -sin(phi)*sin(teta)*ro + cos(phi)*s
		
		
		# velocidade no sistema do CM
		
		vlcmx = gama2*self.vf*ql
		vlcmy = gama2*self.vf*ul
		vlcmz = gama2*self.vf*sl
		
		self.vx = vlcmx+gama1*self.vfx+gama2*gas_atom.vx
		self.vy = vlcmx+gama2*gas_atom.vy
		self.vz = vlcmz+gama1*self.vfz+gama2*gas_atom.vz
				
		self.energy = 0.5*self.mass*self.vf*self.vf
				
		return
		
def ionTrip3D(subject, step, la, f, eField = (0.0, 0.0, 0.0)):
	
	if subject.mass == electronMass:
		charge = -1.60217646E19
	else:
		charge = 1.60217646E19
		
	acceleration = (charge*eField[0]/subject.mass,charge*eField[1]/subject.mass, charge*eField[2]/subject.mass)	

	lastpositionx  = subject.x  
	lastpositiony  = subject.y
	lastpositionz  = subject.z
	
	while True:
		
		subject.vx += acceleration[0]*step
		subject.vy += acceleration[1]*step  
		subject.vz += acceleration[2]*step 
		 
		subject.x = subject.vx*step + subject.x 
		subject.y = subject.vy*step + subject.y
		subject.z = subject.vz*step + subject.z

		distpercorrida = sqrt((subject.x-lastpositionx)**2 + (subject.y-lastpositiony)**2 + (subject.z-lastpositionz)**2)

		# Collision condition: when distance traveled converges to la, collision probability converges to 1.
		if subject.particle_type == "Electron":
			
			if random() < (distpercorrida/la(subject.energy))*0.5 :
				subject.colides()
				lastpositionx = subject.x  
				lastpositiony = subject.y
				lastpositionz  = subject.z
				
			if subject.energy < 1E3:
				f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (subject.x, subject.y,  subject.z, sqrt((subject.x)**2 + (subject.y)**2 + (subject.z)**2), subject.energy, subject.collisioncounter))
				return

		elif subject.particle_type == "Argon+":

			if random() < (distpercorrida/la)*0.5 :
				subject.colides() 	# Change collisions conditions here
				lastpositionx = subject.x  
				lastpositiony = subject.y
				lastpositionz = subject.z
				
			if subject.energy < 1E3:
				f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (subject.x, subject.y, subject.z, sqrt((subject.x)**2 + (subject.y)**2 + (subject.z)**2), subject.energy, subject.collisioncounter))
				return
		else:
			print "panic"
			
			
def simulate_collisions3D(step, free_mean_path,electric_field, ions, particle_type):				
	
	# Some Default values
	energy = 1E9
	x_initial = 0
	y_initial = 0
	z_initial = 0
	
	# Velocity should be normalized to one
	vx_initial = 1
	vy_initial = 0
	vz_initial = 1
	
	f = file('collisionslog3D.txt','w')
	
	f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("#","x", "y", "z", "distance", "energy", "collision"))

	i = 0	
	while i < ions:
		print i
		subject = Particle(particle_type, energy, x_initial, y_initial, z_initial, vx_initial, vy_initial, vz_initial)
		ionTrip3D(subject, step, free_mean_path, f, electric_field) # add eletric field here. default set to 0
		i += 1
	
	f.close()
	
	pass
	
	
if __name__ == '__main__':
		
	# The options are "Argon+" and "Electron"
	particle_type = "Electron"
	
	# Electric field module for x and y directions
	electric_field = (0,0,0) 

	# Measure of effective section
	cross_section = []
	f = open('results.txt', 'r')
	
	for line in f:
		cross_section.append(line[:-1].split('\t'))

	cross_section = map(lambda x: (float(x[0]), float(x[1])),cross_section)
	
	if particle_type == "Electron":

		def free_mean_path(energy):
		
			if energy < cross_section[-1][0]:
				for pair in cross_section:
					if (pair[0]-energy) < 1E-2:
						sigma = pair[1]
			else:
				a = 1.25E-14
				b = 6.72E-1
				sigma = a/(energy**(-b))
			
			return 0.5/(sigma*n)

	else:
		
		free_mean_path = 2.36E-9 # m
			
		
	# Must be very small to be reasonable, according to the free mean path
	step = 1E-15
	
	# Number of ions to be launched
	ions = 1000
	
	print "3D Simulation"
	import profile
	profile.run('simulate_collisions3D(step, free_mean_path, electric_field, ions, particle_type)')
