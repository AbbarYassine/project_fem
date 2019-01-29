import numpy as np
import math
import routinesfem as fem
from ReadMshData import GmshMesh
from createParaview import createVtk
import sys 

# read msh file:

Meshfile = sys.argv[1]
#print(sys.argv[1])

Mesh = GmshMesh(Meshfile)
nodes = Mesh.nodes
elements = Mesh.elements

# functions & params

alpha = (math.pi)

k = 2*math.pi

g_neumann = 0

quad_degre = 2

f = 0

cst = complex(0,1)*k

signe_stiffness = -1.

bord_out_tag = int(sys.argv[4])

bord_in_tag = int(sys.argv[3])

interior = int(sys.argv[2])

def uinc(x,y):
	return -1*np.exp(complex(0,1)*k*(x*np.cos(alpha)+y*np.cos(alpha)))


M = fem.assem_mass(nodes,elements,k,interior)	
D = fem.assem_stiffness(nodes,elements,signe_stiffness,interior)
M_bord = fem.assem_mass_bord(nodes,elements,cst,bord_out_tag)
B = fem.assem_vectorB_dirichlet(nodes,elements,uinc,bord_in_tag)
A = fem.boundary_dirichlet_A(elements,M,M_bord,D,bord_in_tag)
# assembly Matrix 

#A, B = fem.assemble(nodes,elements,f,g_neumann,uinc,k,signe_mass,quad_degre,cst,bord_out_tag,bord_in_tag,interior)
#A = fem.boundary_dirichlet_A(elements,A,bord_in_tag)

 #solve the problem 
u_sol = fem.solve(A,B)

# select triangle from elements
triangles = {k: v for k, v in elements.iteritems() if v[0]==2}
print(A.shape)
print(len(u_sol))
print(len(B))
print(len(triangles))
createVtk(u_sol, nodes, triangles)