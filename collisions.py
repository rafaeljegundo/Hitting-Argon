#!/usr/bin/env python
# encoding: utf-8
"""
collisions2D.py

Created by Rafael Jegundo on 2011-05-11.
Copyright (c) 2011 SMMC. All rights reserved.

Notes: the step needs to be very small, due to the small free mean path. Final results for traveled distance should be ~1E-6

Enunciado: 

Energia cinética inicial ~1MeV. Iões lançados forward.

Depois de cada colisão há isotropia (temos de considerar conservação de momento linear e de energia –
recuo do átomo com que colide, que inicialmente está em repouso)

Densidade do meio – N átomos /cm3 (number density à pressão atmosférica)

A simulação para quando os iões têm energia final de ~1keV. Regista-se a posição final do ião.

Iões de Ar+ em atmosfera de Ar.

Usar como secção eficaz o modelo da esfera rígida (para determinar o caminho livre médio do ião)

Upgrades ao problema (pontos a adicionar resolvida a primeira aproximação):

1) Os átomos do meio têm distribuição de velocidade normal com média em 3/2kT. - DONE
2) Há um campo eléctrico uniforme aplicado.
3) Electrões em vez de iões a serem lançados e uma secção eficaz mais realista que terá de ser
tratada pelo métodos numéricos aprendidos.
4) Passar para 3D.

Test and END.
"""

from random import random
from numpy import pi, sin, cos, sqrt, hypot, arctan
import sys
from time import time


# Constantes Universais

c = 299792458
K = 8.61E-5 # eV/K
T = 293
q = 1.60217646E19 # Coulomb

class Particle:	
	
	"""Ion attributes: position (x,y), velocity (vx,vy) and energy"""
	
	def __init__(self, particle_type, energy, x, y, direction):
		
		argonMass = 931.46E6/(c**2)
		electronMass = 9.1093821E-31
		
		if particle_type == "Argon+":		
			self.mass = argonMass
		elif particle_type == "Electron":
			print "electron!"
			self.mass = electronMass
		
		self.particle_type = particle_type
		self.x = x
		self.y = y
		self.energy = energy
		self.v = sqrt((2*energy)/self.mass)
		self.vx = self.v*cos(direction)
		self.vy = self.v*sin(direction)
				
		self.collisioncounter = 0	
		self.state = ""
		
		return

	def colides(self, collision_type="2D-inerte gas"):

		self.collisioncounter += 1
		self.vi1 = hypot(self.vx,self.vy)
		
		if collision_type == "2D-inerte gas":
			# Válido para electrões ????????????????????????
			teta = 2*pi*random()
			teta1 = arctan(sin(teta)/(1+cos(teta)))
			teta2 = 0.5*(pi-teta)

			self.vf1 = self.vi1/sqrt(1+((sin(teta1)**2)/(sin(teta2)**2)))

		elif collision_type == "2D-3/2KTEnergy gas":
			
			teta = 2*pi*random()
			teta1 = arctan(sin(teta)/(1+cos(teta)))
			teta2 = 0.5*(pi-teta)
			
			direction = random()*180

			gas = Particle("Argon+", (3/2)*K*T, self.x, self.y, direction)
			
			vi2 = hypot(gas.vx, gas.vy)		

			a = 1 + (sin(teta1)/sin(teta2))**2 
			b = 2*vi2*sin(teta1)/sin(teta2)**2*sin(teta) # Será teta? 
			c = -self.vi1**2 + vi2**2 *(-1 + ( sin(teta)/sin(teta2)))
			
			self.vf1 = (-b + sqrt(b*b-4*a*c))/(2*a)
			
		elif collision_type == "3D-3/2KTEnergy gas":
			
			print "TBI"
			
		else:
			
			print "unknown collision type"
			sys.exit()
		
		self.vx = self.vf1*cos(teta1)
		self.vy = self.vf1*sin(teta1)
		self.energy = 0.5*self.mass*self.vf1*self.vf1

		return
		
	def log(self):
		self.state = map(str,[self.x, self.y, self.energy, self.vx, self.vy, self.collisioncounter])
		return

def ionTrip2D(subject, step, la, f, eField = (0.0,0.0)):
	
	acceleration = (q*eField[0]/subject.mass,q*eField[1]/subject.mass)	
	
	while True:
		lastpositionx  = subject.x  
		lastpositiony  = subject.y
		subject.vx += acceleration[0]*step
		subject.vy += acceleration[1]*step  
		subject.x = subject.vx*step + subject.x 
		subject.y = subject.vy*step + subject.y
		distpercorrida = sqrt((subject.x-lastpositionx)**2 + (subject.y-lastpositiony)**2)
		
		# Collision condition: when distance traveled converges to la, collision probability converges to 1.
		# Deverá evoluir para colisão nula
		
		if random() < distpercorrida/la*0.5 :
			subject.colides()	
			subject.log()
		if subject.energy < 1E3:
			f.write("%s\t%s\t%s\t%s\t%s\n" % (subject.x, subject.y, hypot(subject.x,subject.y), subject.energy, subject.collisioncounter))
			return
			
	
def simulate_collisions2D(step, free_mean_path,electric_field, ions, particle_type):

	# Some Default values
	energy = 1E9
	x_initial = 0
	y_initial = 0
	direction = 0 # Degrees between ion direction and the positive X axe.
	
	f = file('collisionslog2D.txt','w')
	
	f.write("%s\t%s\t%s\t%s\t%s\n" % ("x", "y", "distance", "energy", "collision"))

	i = 0	
	while i < ions:
		subject = Particle(particle_type, energy, x_initial, y_initial, direction)
		ionTrip2D(subject, step, free_mean_path, f, electric_field) # add eletric field here. default set to 0
		i += 1
	
	f.close()
	
	pass
	



if __name__ == '__main__':
	print u"Falta a probabilidade de colisão nula"

	print u"Note-se que no caso de electrões, e = -q"
	
	# The options are "Argon+" and "Electron"
	particle_type = "Argon+"
	
	# Electric field module for x and y directions
	electric_field = (0,0) 

	# Measure of effective section
	free_mean_path = 2.36E-9 
	
	# Must be very small to be reasonable, according to the free mean path
	step = 0.00000000001
	
	# Number of ions to be launched
	ions = 20

	print "2D Simulation"
	simulate_collisions2D(step, free_mean_path, electric_field, ions, particle_type)
	