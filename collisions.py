#!/usr/bin/env python
# encoding: utf-8
"""
main.py

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

1) Os átomos do meio têm distribuição de velocidade normal com média em 3/2kT.
2) Há um campo eléctrico uniforme aplicado.
3) Electrões em vez de iões a serem lançados e uma secção eficaz mais realista que terá de ser
tratada pelo métodos numéricos aprendidos????????.
4) Passar para 3D.

Test and END.
"""

from random import random
from numpy import pi, sin, cos, sqrt, hypot
import math
from time import time


# Constantes Universais

c = 299792458
K = 8.61E-5 # eV/K
T = 293
q = 1.60217646E19 # Coulomb

class Particle:	
	
	"""Ion attributes: position (x,y), velocity (vx,vy) and energy"""
	
	def __init__(self, tipe, x, y, vx, vy, energy):
		
		argonMass = 931.46E6/(c**2)
		electronMass = 9.1093821E-31
		
		self.tipe = tipe
		self.x = x
		self.y = y
		self.vx = vx
		self.vy = vy		
		self.energy = energy
		self.collisioncounter = 0	
		self.state = ""
		
		if tipe == "Argon+":		
			self.mass = argonMass
		elif self.type == "Electron":
			print "electron!"
			self.mass = electronMass
		return
		
	def colides(self):
		
		velo = sqrt(((3/2)*K*T)/(2*self.mass))
		vxrandom = sqrt(random()*velo)
		vyrandom = sqrt(velo**2-vxrandom**2)
		
		gas = Particle(self.tipe, self.x,self.y,vxrandom,vyrandom, (3/2)*K*T)

		# Falta fazer as contas
		
		self.collisioncounter += 1
		
		teta = random()*2*pi
		teta1 = math.atan(sin(teta)/(1+cos(teta)))
		teta2 = 0.5*(pi-teta)
		
		# velocidade
		vi2 = hypot(gas.vx, gas.vy)
		self.vi1 = hypot(self.vx,self.vy)
		
		a = 1 - (sin(teta1)*sin(teta1))/(sin(teta2)*sin(teta2))
		b = 2 * (vi2*sin(teta1)*sin(teta))/(sin(teta2)*sin(teta2))
		c = ((vi2**2)*sin(teta)**2)/sin(teta2)
				
		try:
			self.vf1 = (-b + math.sqrt(b*b-4*a*c))/(2*a)
		except:
			try:
				self.vf1 = (-b - math.sqrt(b*b-4*a*c))/(2*a)
			except:
				self.vf1 = self.vi1*(1/sqrt(1+(sin(teta1)**2)/(sin(teta2)**2))) # Original
				print "used fail-safe. It shouldn't work this way :\\ To Be corrected"
			
		self.vx = self.vf1*cos(teta1)
		self.vy = self.vf1*sin(teta1)

	    #energia
		if self.tipe == "Argon+":
			self.energy = self.energy*0.5 # 0.5 é o valor espectável relativo
		else:
			self.energy = 0.5*self.mass*self.vf1*self.vf1
		
		return
		
	def log(self):
		self.state = map(str,[self.x, self.y, self.energy, self.vx, self.vy, self.collisioncounter])
		return

def ionTrip(subject, step, la, f, eField = (0.0,0.0)):
	acceleration = (q*eField[0]/subject.mass,q*eField[1]/subject.mass)	
	while True:
		lastpositionx  = subject.x  
		lastpositiony  = subject.y
		subject.vx += acceleration[0]*step
		subject.vy += acceleration[1]*step  
		subject.x = subject.vx*step + subject.x 
		subject.y = subject.vy*step + subject.y
		distpercorrida = sqrt((subject.x-lastpositionx)**2 + (subject.y-lastpositiony)**2)
		if random() < (distpercorrida/la)*0.5 :
			subject.colides()	
			subject.log()
			print "colision"
		if subject.energy < 1E3:
			f.write("%s\t%s\t%s\t%s\t%s\n" % (subject.x, subject.y, hypot(subject.x,subject.y), subject.energy, subject.collisioncounter))
			return
	
	
def main():
	
	ions = 10
	i = 0
	step = 0.000000001 # tem de ser muito pequeno para ter validade neste contexto (lambda = la = 2E-9)
	la = 2.36E-9
#	eField = (10,10) # Electric field module for x and y directions
	
	f = file('collisionslog.txt','w')
	f.write("%s\t%s\t%s\t%s\t%s\n" % ("x", "y", "distance", "energy", "collision"))
	
	while i < ions:
		subject = Particle("Argon+",0,0,100,3,1E9)
		print "new subject"
		ionTrip(subject, step, la, f) # add eletric field here. default set to 0
		i += 1
	f.close()
	pass


if __name__ == '__main__':
	main()

