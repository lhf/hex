# A compact representation for adaptive hexagonal meshes
#
# usage: python hex.py R INPUT -- all arguments optional

from __future__ import print_function
from sys import argv
from random import random

# depth of refinement
R= int(argv[1]) if len(argv)>1 else 6

# file containing initial mesh
INPUT= argv[2] if len(argv)>2 else "-"

# size of base mesh
N=5
L=1

# make array
def array(n):
	return list(range(n))

# lattice coordinates for basis 1+w2, -1+2*w2
def lattice(a,b):
	return complex(a,b)

# cartesian coordinates for basis 1+w2, -1+2*w2
def cartesian(z):
	w2=complex(0.5,0.86602540378443864676372317075293618347140262690519031)
	a,b=z.real,z.imag
	return L*((a-b)+(a+2*b)*w2)	# = L*(a*(1+w2)+b*(-1+2*w2))

# scale and translate basic anchors
def A(c,s,k):
	return c+STAR[k]/(2**s)

# face vertex
def vertex(c,k):
	t=F[c].t
	s=F[c].s
	if t==6:
		return cartesian(c)+VERT[k]/(2**s)
	else:
		return cartesian(c)+QUAD[t][k]/(2**s)-cartesian(A(0,s+1,t))

# face adjacent to quad along longest edge
def mate(c):
	t=F[c].t
	s=F[c].s
	return A(c,s+1,t)

# basic star of hex centers
STAR=array(6+1)
STAR[0]=lattice(1,0)
STAR[1]=lattice(0,1)
STAR[2]=lattice(-1,1)
STAR[3]=lattice(-1,0)
STAR[4]=lattice(0,-1)
STAR[5]=lattice(1,-1)
STAR[6]=STAR[0]

# basic hex vertices
VERT=array(6)
for k in range(6):
	VERT[k]=cartesian((0+STAR[k]+STAR[k+1])/3)

# basic quad vertices
QUAD=array(6)
for k in range(6):
	QUAD[k]=[VERT[k],VERT[k]/2,VERT[k-1]/2,VERT[k-1]]

# use dot notation for dict -- https://stackoverflow.com/a/74214556/107090
class DOTTED(dict):
	__getattr__ = dict.get
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__

# faces
F={}
Q=set()

def addface(c,t,s):
	if c in F:
		F[c].t=6
		F[c].s=F[c].s+1
	else:
		F[c]=DOTTED({'t':t,'s':s})
	Q.add(c)

def baseface(c):
	addface(c,6,0)

# base mesh: hexagonal grid
def basemesh(N):
	b=lattice(0,0)
	for i in range(N):
		for j in range(N):
			c=b+lattice(0,j)
			baseface(c)
			c=c+lattice(1,0)
			baseface(c)
		b=b+lattice(2,-1)

# load mesh from csv file
def loadmesh(filename):
	N=0
	for line in open(filename):
		if N>0:
			a,b,t,s=list(map(float,line.split(",")))
			c=complex(a,b)
			t=int(t)
			s=int(s)
			addface(c,t,s)
		N=N+1

# subdivide hex
def subdivide(c):
	s=F[c].s
	for k in range(6):
		h=A(c,s+1,k)
		addface(h,k,s)
	F[c].s=s+1

# subdivide quad at border
def subdivide4(c):
	t=F[c].t
	s=F[c].s
	for j in range(2,4+1):
		k=(t+j)%6
		h=A(c,s+2,k)
		addface(h,k,s+1)
	F[c].s=s+1

# implicit curve (Taubin, 1994)
def f(v):
	z=v-complex(7,3)
	z=0.5*z
	x,y=z.real,z.imag
	z=0.004+0.110*x-0.177*y-0.174*x*x+0.224*x*y-0.303*y*y-0.168*x*x*x+0.327*x*x*y-0.087*x*y*y-0.013*y*y*y+0.235*x*x*x*x-0.667*x*x*x*y+0.745*x*x*y*y-0.029*x*y*y*y+0.072*y*y*y*y;
	return z

# refinement criterion: uniform
def needsrefinement(v):
	return F[v].t==6 and F[v].s < R

# refinement criterion: random
def needsrefinement(v):
	return F[v].t==6 and F[v].s < R and random() < 0.55

# refinement criterion: implicit curve
def needsrefinement(v):
	n = 6 if F[v].t==6 else 4
	w = [f(vertex(v,k)) for k in range(n)]
	w.append(w[0])
	z = min([w[k]*w[k+1] for k in range(n)])
	return F[v].s < R and z<=0

# refine respecting border
def refine(c):
	if F[c].t==6:
		subdivide(c)
		Q.add(c)
	else:
		m=mate(c)
		if not (m in F):
			subdivide4(c)
		else:
			refine(m)
			assert(F[c].t==6)
			refine(c)

# refine extending border
# does not work on unbounded implicit curves
def refine(c):
	if F[c].t==6:
		subdivide(c)
		Q.add(c)
	else:
		m=mate(c)
		if not (m in F):
			addface(m,6,F[c].s)
		refine(m)
		assert(F[c].t==6)
		refine(c)

def drawface(c):
	print("")
	print("% face centered at",c,F[c])
	n = 6 if F[c].t==6 else 4
	w="moveto"
	for k in range(n):
		z=vertex(c,k)
		x,y=z.real,z.imag
		print(x,y,w)
		w="lineto"
	if needsrefinement(c):
		print("b"+str(n))
	else:
		print("f"+str(n))
	z=cartesian(c)
	x,y=z.real,z.imag
	print("%",x,y,"c"+str(n))

def dual(c,s,k):
	d=A(c,s+1,k)
	if d in F:
		return d
	d=A(c,s,k)
	if d in F:
		return d
	else:
		return None

def dualedge(c,s,k):
	k=k%6
	d=dual(c,s,k)
	#print "% edge",k,c,d
	if d!=None:
		z1=cartesian(c)
		z2=cartesian(d)
		t=0.0; z3=(1-t)*z1+t*z2; x3,y3=z3.real,z3.imag
		t=1-t; z4=(1-t)*z1+t*z2; x4,y4=z4.real,z4.imag
		if x3<x4 or (x3==x4 and y3<y4):
			print(x3,y3,x4,y4,"moveto lineto stroke")

def drawdual(c):
	s=F[c].s
	t=F[c].t
	if t==6:
		#continue
		for k in range(6):
			dualedge(c,s,k)
	else:
		dualedge(c,s,t)
		dualedge(c,s+1,t+2)
		dualedge(c,s+1,t+3)
		dualedge(c,s+1,t+4)

# main -----------------------------------------------------------------------

# initial mesh
if INPUT=="-":
	basemesh(N)
else:
	loadmesh(INPUT)

# refine mesh
while Q:
	c=next(iter(Q))
	Q.remove(c)
	if needsrefinement(c):
		refine(c)

# output PostScript
print('''\
%!PS-Adobe-2.0 EPSF-2.0
%%BoundingBox: 18 -4 487 337
50 50 translate
30 dup scale
0.05 setlinewidth
0 setlinewidth
1 setlinejoin
1 setlinecap
/c { 0.05 0 360 arc fill } bind def
/v { 0.1 0 360 arc fill } bind def
/f { closepath gsave fill grestore 0 0 0 setrgbcolor stroke } bind def
/f4 { 0.8 0.8 1 setrgbcolor f } bind def
/f6 { 0.8 1 0.8 setrgbcolor f } bind def
/c4 { 0 0 1 setrgbcolor c } bind def
/c6 { 0 1 0 setrgbcolor c } bind def
/b4 { 1 0 0 setrgbcolor f } bind def
/b6 { 1 0 0 setrgbcolor f } bind def
/f4 { 1 0.8 0.5 setrgbcolor f } bind def
/f6 { 1 0.8 0.5 setrgbcolor f } bind def
%/b4 { 1 0.8 0.5 setrgbcolor f } bind def
%/b6 { 1 0.8 0.5 setrgbcolor f } bind def
''')
print("%=ARG",R,INPUT)

R+=2	# identify leaf faces
print("")
print("% faces")
for c in F:
	drawface(c)
	pass

print("")
print("% dual")
print("0 0.6 0 setrgbcolor")
print("0.01 setlinewidth")
for c in F:
	#drawdual(c)
	pass

print("")
print("showpage")
print("%%EOF")

#exit()

# output CSV
print("")
PREFIX="%=CSV"
print(PREFIX, "a,b,t,s")
for c in F:
	a,b=c.real,c.imag
	print(PREFIX, ','.join(map(str,[a,b,F[c].t,F[c].s])))

