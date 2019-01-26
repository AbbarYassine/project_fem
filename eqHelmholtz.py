import numpy as np
import math
import routinesfem as fem
from ReadMshData import GmshMesh

# read msh file:

Mesh = GmshMesh('disques.msh')
nodes = Mesh.nodes
elements = Mesh.elements

# functions & params

alpha = math.pi

k = 2*math.pi

g_neumann = 0

quad_degre = 2

f = 0

cst = complex(0,1)*k

signe_mass = -1.

bound_out = 12

bound_in = 11

interior = 10

def uinc(x,y):
	return -1*np.exp(complex(0,1)*k*(x*np.cos(alpha)+y*np.cos(alpha)))

# assembly Matrix :	

A, B = fem.assemble(nodes,elements,f,g_neumann,k,signe_mass,quad_degre,cst,bound_out,bound_in,interior)
u_sol = fem.solve(A,B)
print(u_sol)

	
