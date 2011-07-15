#!/usr/bin/env python
# encoding: utf-8

from random import random
from numpy import pi, sin, cos, sqrt, hypot, arctan
from sys import exit
from time import time


# Constantes Universais

c = 299792458 # m/s
K = 8.61E-5 # eV/K
T = 293 # K
argonMass = 931.46E6/(c**2) # eV
electronMass = 511E3/(c**2) # eV
n = 2.18E19 # àtomos /cm^3


class Particle:	
	
	"""Ion attributes: position (x,y), velocity (vx,vy) and energy"""
	
	def __init__(self, particle_type, energy, x, y, direction):
		
		self.collisioncounter = 0	
		
		if particle_type == "Argon+":		
			self.mass = argonMass
		elif particle_type == "Electron":
			self.mass = electronMass
		
		self.particle_type = particle_type

		self.x = x
		self.y = y

		self.energy = energy

		self.v = sqrt((2*energy)/self.mass)
		self.vx = self.v*cos(direction)
		self.vy = self.v*sin(direction)
				
		return

	def colides(self, collision_type="3/2KT Energy gas", energy_method="f(v)"):

		self.collisioncounter += 1
		
		self.vi1 = hypot(self.vx,self.vy)
		
		if self.particle_type == "Argon+":
			
			teta = 2*pi*random()
			teta1 = arctan(sin(teta)/(1+cos(teta)))
			teta2 = 0.5*(pi-teta)
			
		elif self.particle_type == "Electron":
			
			teta = 2*pi*random()
			teta1 = arctan(sin(teta)/((electronMass/argonMass)+cos(teta)))
			teta2 = 0.5*(pi-teta)
		
		# Exit velocity calculation for each collision type
		if collision_type == "0 Energy gas":
 
			self.vf1 = self.vi1/sqrt(1+((sin(teta1)**2)/(sin(teta2)**2)))

		elif collision_type == "3/2KT Energy gas":
			
			direction = random()*pi

			gas = Particle("Argon+", (3/2)*K*T, self.x, self.y, direction)
			
			vi2 = hypot(gas.vx, gas.vy)		

			a = 1 + (sin(teta1)/sin(teta2))**2 
			b = (2*vi2*sin(teta1)*sin(teta))/sin(teta2)**2 
			c = -self.vi1**2 +(vi2**2)*(-1 + (sin(teta)/sin(teta2)))
			
			self.vf1 = (-b + sqrt(b*b-4*a*c))/(2*a)
			
		elif collision_type == "different mass":
			
			direction = random()*pi

			gas = Particle("Argon+", (3/2)*K*T, self.x, self.y, direction)
			
			vi2 = hypot(gas.vx, gas.vy)		

			m2 = self.mass
			m1 = gas.mass
			
			a = m1 + m2*(sin(teta1)/sin(teta2))**2 
			b = (2*vi2*sin(teta1)*sin(teta)*m2)/sin(teta2)**2 
			c = -m1*self.vi1**2 + (vi2**2)*m2*(-1 + ( sin(teta)/sin(teta2)))
			
			self.vf1 = (-b + sqrt(b*b-4*a*c))/(2*a)
			
		else:
			
			print "unknown collision type"
			exit()
		
		# Final velocity components calculations
		self.vx = self.vf1*cos(teta1)
		self.vy = self.vf1*sin(teta1)
				
		# Final energy calculations
		if energy_method == "f(v)":
			self.energy = 0.5*self.mass*self.vf1*self.vf1
 		elif energy_method == "half":
			self.energy = 0.5*self.energy
		else:
			print "unknown energy calculation method"
			exit()
		return

def ionTrip2D(subject, step, la, f, eField = (0.0,0.0)):
		
	if subject.mass == electronMass:
		charge = -1.60217646E19
	else:
		charge = 1.60217646E19
	
	acceleration = (charge*eField[0]/subject.mass,charge*eField[1]/subject.mass)
	
	lastpositionx  = subject.x  
	lastpositiony  = subject.y
	
	while True:
		subject.vx += acceleration[0]*step
		subject.vy += acceleration[1]*step  
		subject.x += subject.vx*step
		subject.y += subject.vy*step
		distpercorrida = sqrt((subject.x-lastpositionx)**2 + (subject.y-lastpositiony)**2)
		# Collision condition: when distance traveled converges to la, collision probability converges to 1.
		
		if subject.particle_type == "Electron":
			
			if random() < (distpercorrida/la(subject.energy))*0.5 :
				subject.colides("different mass", "f(v)")	# Change collisions conditions here
				lastpositionx = subject.x  
				lastpositiony = subject.y
				
			if subject.energy < 1E3 :
				f.write("%s\t%s\t%s\t%s\t%s\n" % (subject.x, subject.y, hypot(subject.x,subject.y), subject.energy, subject.collisioncounter))
				return

		elif subject.particle_type == "Argon+":

			if random() < (distpercorrida/la)*0.5 :
				subject.colides("3/2KT Energy gas", "f(v)") 	# Change collisions conditions here
				lastpositionx = subject.x  
				lastpositiony = subject.y
				
			if subject.energy < 1E3 :
				f.write("%s\t%s\t%s\t%s\t%s\n" % (subject.x, subject.y, hypot(subject.x,subject.y), subject.energy, subject.collisioncounter))
				return
		else:
			print "panic"
	return
			
	
def simulate_collisions2D(step, free_mean_path,electric_field, ions, particle_type):

	# Some Default values
	energy = 1E6 #eV
	x_initial = 0
	y_initial = 0
	direction = 0 # Degrees between ion/electron direction and the positive X axe.
	
	f = file('collisionslog2D.txt','w')
	
	f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("#", "x", "y", "distance", "energy", "collision"))

	i = 0	
	while i < ions:
		print i
		subject = Particle(particle_type, energy, x_initial, y_initial, direction)
		ionTrip2D(subject, step, free_mean_path, f, electric_field) # add eletric field here. default set to 0
		i += 1
	
	f.close()
	
	pass
	

if __name__ == '__main__':
	
	# The options are "Argon+" and "Electron"
	particle_type = "Argon+"
	
	# Electric field module for x and y directions
	electric_field = (-0.001,0) 

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
		
	# Must be << free_mean_path
	step = 1E-15
	
	# Number of ions to be launched
	ions = 1000
			
	print "2D Simulation"
	import profile
	profile.run('simulate_collisions2D(step, free_mean_path, electric_field, ions, particle_type)')
	