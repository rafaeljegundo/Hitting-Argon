#!/usr/bin/env python
# encoding: utf-8
"""
main.py

Created by Rafael Jegundo on 2011-05-11.
Copyright (c) 2011 SMMC. All rights reserved.
"""

from random import random
from math import pi, atan, sin, cos, sqrt, hypot
from time import time

# Constantes Universais

c = 299792458


class Ion:	
	
	"""Ion attributes: position (x,y), velocity (vx,vy) and energy"""
	
	def __init__(self,x,y,vx,vy,energy):
		
		self.x = x
		self.y = y
		self.vx = vx
		self.vy = vy
		self.mass = 931.46E6/(c**2) # useless for now
		self.energy = energy
		self.collisioncounter = 0	
		self.state = ""
		return
		
	def colides(self):
		
		self.collisioncounter += 1
		
		teta = random()*2*pi
		teta1 = atan(sin(teta)/(1+cos(teta)))
		teta2 = 0.5*(pi-teta)
		
		# velocidade
		self.vi = hypot(self.vx,self.vy)
		self.vf = self.vi*(1/sqrt(1+(sin(teta1)**2)/(sin(teta2)**2)))
		self.vx = self.vf*cos(teta1)
		self.vy = self.vf*sin(teta1)

	    # energia
		self.energy = self.energy*0.5 # 0.5 é o valor espectável relativo
		return
		
	def log(self):
		self.state = map(str,[self.x, self.y, self.energy, self.vx, self.vy, self.collisioncounter])
		return

def ionTrip(subject, step, la, f):
	while True:
		lastpositionx  = subject.x  
		lastpositiony  = subject.y  
		subject.x  = subject.vx*step  +  subject.x 
		subject.y  =  subject.vy*step +  subject.y
		distpercorrida = sqrt((subject.x-lastpositionx)**2 + (subject.y-lastpositiony)**2)
		if random() < (distpercorrida/la)*0.5 :
			subject.colides()	
			subject.log()
		if subject.energy < 1E3:
			f.write("%s\t%s\t%s\n" % (hypot(subject.x,subject.y),subject.energy,subject.collisioncounter))
			return
	
	
def main():
	
	ions = 10
	i = 0
	step = 0.00000001 # tem de ser muito pequeno para ter validade neste contexto (lambda = la = 2E-9)
	la = 2.36E-9
	
	# Ficheiro: distance	energy	collisions
	f = file('collisionslog.txt','w')
	
	t1 = time()
	while i < ions:
		subject = Ion(0,0,100,3,1E9)
		ionTrip(subject, step, la, f)
		i += 1
	f.close()
	t2 = time()
	print t2-t1
	pass


if __name__ == '__main__':
	main()

