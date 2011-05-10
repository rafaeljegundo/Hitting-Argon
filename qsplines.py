#!/usr/bin/env python
# encoding: utf-8
"""
qsplines.py

Uses tridiagonal.py from http://www.math-cs.gordon.edu/courses/mat342/python.html

Created by Rafael Jegundo on 2011-05-10.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import tridiagonal

def main():
	
	f = open('data_acquired.txt','r')
	data = []
	x = []
	y = []
	h = []
	u = []
	v = []
	
	
	for line in f:
		data.append(map(eval,line.split('\t')))

	n = len(data) - 2
	
	for pair in data:
	  x.append(pair[0])
	  y.append(pair[1])

	
	i = 0
	while i < n:
		h.append(x[i+1]-x[i])
		i += 1
		
	i = 0
	while i < n: 
		try:
			u.append(2*(h[i] + h[i-1]))
		except IndexError:
			u.append(0) # avaliar depois com os primeiros resultados. 0 Vs h[i]
		i+=1

	i = 0 
	while i < n: 
		v.append(((6/h[i])*(y[i+1]-y[i]))-(6/(h[i-1]))*(y[i]-y[i-1]))
		i += 1
		
		
	
	r = open('results.txt','w')
	
	pass


if __name__ == '__main__':
	main()

