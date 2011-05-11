#!/usr/bin/env python
# encoding: utf-8
"""
qsplines.py

Created by Rafael Jegundo on 2011-05-10.
Copyright (c) 2011 SMMC. All rights reserved.
"""
import numpy as np
import matplotlib
matplotlib.use('macosx')

import pylab as plt

def tridiagonal(a,b,c,v):
	
	n = len(b)
	
	for i in range(1,n-1):
	    m = a[i]/b[i-1]
	    b[i] -= m*c[i-1]
	    v[i] -= m*v[i-1]

	x = np.zeros(n)
	x[-1] = v[-1]/b[-1]

	print x,n
	print len(x), len(v), len(c), len(b)
	
	for i in range(n-2,-1,-1):
	    x[i] = (v[i]-c[i]*x[i+1])/b[i]
		
	return x

def splinefun(x,y,f,h,xr):
	
	yr = []
	n = len(f)
	pontos = len(xr)
	print len(h), len(f), len(x), len(y),n
	print len(xr), pontos
	
	for p in range(pontos):
			for a in range(n-1):
				if x[a] <= xr[p] and x[a+1] >= xr[p]:
					print a,p, x[a], xr[p], x[a+1]
					yr.append( \
						(f[a]/(6*h[a]))*(x[a+1]-xr[p])**3 +\
						(f[a+1]/(6*h[a]))*(xr[p]-x[a])**3	+\
						((y[a+1]/h[a])-((f[a+1]*h[a])/6))*(xr[p]-x[a]) +\
						((y[a]/h[a])-(((f[a]*h[a])/6))*(x[a+1]-xr[a])))
				else:
					continue
			
	return yr
	
	
def main():
	
	# Gathering original data
	f = open('data_acquired.txt','r')
	data = []
	for line in f:
		data.append(map(eval,line.split('\t')))
	
	x = []
	y = []
	u = []
	h = []
	v = []
	n = len(data) - 2
	
	# Defining x and y vectors from original data
	for pair in data:
	  x.append(pair[0])
	  y.append(pair[1])

	# Calcular u, h e v
	
	i = 0
	while i < n:
		h.append(x[i+1]-x[i])
		i += 1
		
	i = 0
	while i < n: 
		try:
			u.append(2*(h[i] + h[i-1]))
		except IndexError:
			u.append(0)
		i+=1

	i = 0 
	while i < n: 
		v.append(((6/h[i])*(y[i+1]-y[i]))-(6/(h[i-1]))*(y[i]-y[i-1]))
		i += 1
	
	xr = np.linspace(1E-4,50,100)	
	
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
		
	# Another hack. The problem is on yr calculus.
	print len(x), len(y)
	print len(xr[:88]), len(yr)
	plt.plot(x,y,'ro',xr[:88],yr,'bo')
	plt.show()
	
	pass


if __name__ == '__main__':
	main()

