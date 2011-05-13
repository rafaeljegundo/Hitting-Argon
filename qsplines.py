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

	
	for i in range(n-2,-1,-1):
	    x[i] = (v[i]-c[i]*x[i+1])/b[i]
		
	return x

def splinefun(x,y,f,h,xr):
	
	yr = []
	n = len(f)
	
	for p in range(len(xr)-1):
			for a in range(n-1):
				if x[a] <= xr[p] and x[a+1] >= xr[p]:
					#ai = (1/6)*((f[a+1]-f[a])/h[a])
					#bi = f[a]/2
					#ci = (f[a+1]-f[a])/h[a] - (h[a]*f[a+1] + 2*h[a]*f[a])/6
					#try:
						#di = yr[-1]
					#except:
						#print "one"
						#di = 0
					c = y[a+1]/h[a] - (h[a]/6)*f[a+1]
					d = y[a]/h[a] - (h[a]/6)*f[a]
					yr.append((f[a]/(6*h[a]))*(x[a+1]-xr[p])**3 +\
						(f[a+1]/(6*h[a]))*(xr[p]-x[a])**3 +\
						c*(xr[p]-x[a]) +\
						d*(x[a+1]-xr[p]))
					print a,p, x[a], xr[p], x[a+1], yr[p], y[a]		
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
	
	n = len(data)
	
	# Defining x and y vectors from original data
	for pair in data:
	  x.append(pair[0])
	  y.append(pair[1])

	# Calcular u, h e v
	
	for i in range(n-1):
		h.append(x[i+1]-x[i])
		
	for i in range(n-1):
		u.append(2*(h[i] + h[i-1]))
	
	for i in range(n-1):
		v.append(((6/h[i])*(y[i+1]-y[i]))-(6/(h[i-1]))*(y[i]-y[i-1]))
	
	xr =[]
	
	# Interpolating one point between each knott
	
	for i in range(len(x)-1):
		xr.append(x[i]+(x[i+1]-x[i])/2)
	
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
	
	plt.plot(x,y,'ro',xr[:-1],yr,'bo') # HAS BUG -1
	plt.show()
	
	pass


if __name__ == '__main__':
	main()

