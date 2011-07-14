#!/usr/bin/env python
# encoding: utf-8
"""
qsplines.py

Created by Rafael Jegundo on 2011-05-10.
Copyright (c) 2011 SMMC. All rights reserved.
"""
import numpy as np
from random import random

def tridiagonal(a,b,c,v):
	
	n = len(b)
	
	for i in range(1,n-1):
	    m = a[i]/b[i-1]
	    b[i] -= m*c[i-1]
	    v[i] -= m*v[i-1]

	x = np.zeros(n)
	x[-1] = v[-1]/b[-1]

	
	for i in range(n-2,-1,-1):
	    x[i] = (v[i]-c[i]*x[i+1])/b[i]
		
	return x

def splinefun(x,y,f,h,xr):
	
	yr = []
	n = len(f)
	p = 0
	
	#print len(f), len(x), len(y), len(h), len(xr)
	
	for p in range(len(xr)):
			for a in range(n-1):
				
				
				if (x[a] <= xr[p]) and (x[a+1] >= xr[p]):
					
					#print p,a,x[a], xr[p], x[a+1]
					
					c = y[a+1]/h[a] - (h[a]/6)*f[a+1]
					
					d = y[a]/h[a] - (h[a]/6)*f[a]
					
					yr.append((f[a]/(6*h[a]))*(x[a+1]-xr[p])**3 +\
						(f[a+1]/(6*h[a]))*(xr[p]-x[a])**3 +\
						c*(xr[p]-x[a]) +\
						d*(x[a+1]-xr[p]))
					
					#print a,p, x[a], xr[p], x[a+1], yr[p], y[a]	
					
				else:
					continue
			p +=1
	return yr
	
	
def interpolspline(x,y):
		
	u = []
	h = []
	v = []
	
	n = len(x)

	# Calcular u, h e v
	
	for i in range(n-1):
		h.append(x[i+1]-x[i])
		
	for i in range(n-1):
		u.append(2*(h[i] + h[i-1]))
	
	for i in range(n-1):
		v.append(((6/h[i])*(y[i+1]-y[i]))-(6/(h[i-1]))*(y[i]-y[i-1]))
	
	xr = np.linspace(x[0],x[-1],10000)
		
	f = tridiagonal(h,u,h,v)
	
	yr = splinefun(x,y,f,h,xr)
	
	r = open('results.txt','w')

	# writing to file	
	i = 0
	while True:
		try:
			r.write(str(xr[i]) + '\t' + str(yr[i]) + '\n')
			i+=1
		except IndexError:
			break
			
	r.close()
	
	pass


if __name__ == '__main__':
	
	# Just a test
	
	f = open('cross_section.txt','r')
	
	data = []
	
	for line in f:
		data.append(map(eval,line.split('\t')))
	
	x = []
	y = []
	
	# Defining x and y vectors from original data
	for pair in data:
	  x.append(pair[0])
	  y.append(pair[1])
	
	interpolspline(x,y)

