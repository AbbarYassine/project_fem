from __future__ import absolute_import, division, print_function
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


def quadrature_triangle():
    """Gauss points for a triangle element (3 points)
    Returns : W weights and Xp points
    """
    W = np.zeros([3])
    Xp = np.zeros([3, 2])
    W[0] = 1./6
    W[1] = 1./6
    W[2] = 1./6

    Xp[0, 0] = 1./6
    Xp[1, 0] = 4./6
    Xp[2, 0] = 1./6

    Xp[0, 1] = 1./6
    Xp[1, 1] = 1./6
    Xp[2, 1] = 4./6

    return W, Xp
    
def quadrature_segment(num_points):
    """Return points and weights for Newton-Cotes rules on [-1,1]."""
    n = num_points
    points  = np.zeros(n)
    weights = np.zeros(n)
    if n == 1:
        # Midpoint rule
        points[0] = 0
        weights[0] = 2
    elif n == 2:
        # Trapezoidal rule
        points[:] = [-1, 1]
        weights[:] = [1, 1]
    elif n == 3:
        # Simpson's rule
        points[:] = [-1, 0, 1]
        weights[:] = [1./3, 4./3, 1/3.]
    elif n == 4:
        # Simpson's 3/8 rule
        points[:] = [-1, -1./3, 1./3, 1]
        weights[:] = [1./4, 3./4, 3/4., 1./4]
    else:
        raise ValueError('Newton-Cotes formula with %d>4 not done' % n)
        
    return points, weights


def elemental_mass_matrix(node1,node2,node3,k):
    """ Return 3x3 elemental mass matrix of a triangle element """
    detJk = abs((node2[0]-node1[0])*(node3[1]-node1[1]) - (node3[0]-node1[0])*(node2[1]-node1[1]))
    Mref = np.matrix('2. 1. 1.; 1. 2. 1.; 1. 1. 2.')
    Me = k**2*detJk*Mref 
    return Me     
 
def elemental_stiffness_matrix(node1,node2,node3,signe):
    """ return 3*3 elemental striffness matrix"""
    detJk = abs((node2[0]-node1[0])*(node3[1]-node1[1]) - (node3[0]-node1[0])*(node2[1]-node1[1]))
    Kp = 0.5*detJk
    deltaphiref = np.matrix('-1. 1. 0.0;-1. 0. 1.')
    a = node3[1]-node1[1]
    b = node1[1]-node2[1]
    c = node1[0]-node3[0]
    d = node2[0]-node1[0]
    Bkp = 1./detJk*np.matrix([[a,b],[c,d]])
    Dep = np.matrix('0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0')
    for i in range(3):
        for j in range(3):
            Dep[i,j] = signe*Kp*(deltaphiref[:,j].T)*(Bkp.T)*(Bkp.H.T)*(deltaphiref[:,i].H.T)
    return Dep        

	    

def elemental_vector(node1,node2,node3,f):
    """ Return 3x1 f elemental vector of a triangle element """
    detJk = abs((node2[0]-node1[0])*(node3[1]-node1[1]) - (node3[0]-node1[0])*(node2[1]-node1[1]))
    bp = np.matrix('0.0 0.0 0.0').T
    node1 = node1[0:2]
    node2 = node2[0:2]
    node3 = node3[0:2]
    W, Xp = quadrature_triangle()

    for i in range(3):
        xm = (1-Xp[i,0] -Xp[i,1])*np.array(node1) + Xp[i,0]*np.array(node2) + Xp[i,1]*np.array(node3)
        bp [0] +=W[i]*f(xm[0],xm[1])*(1-Xp[i,0] -Xp[i,1])

    for i in range(3):
        xm = (1-Xp[i,0] -Xp[i,1])*np.array(node1) + Xp[i,0]*np.array(node2) + Xp[i,1]*np.array(node3)
        bp [1] +=W[i]*f(xm[0],xm[1])*(Xp[i,0])

    for i in range(3):
        xm = (1-Xp[i,0] -Xp[i,1])*np.array(node1) + Xp[i,0]*np.array(node2) + Xp[i,1]*np.array(node3)
        bp [2] +=W[i]*f(xm[0],xm[1])*(Xp[i,1])
        
    bp = detJk*bp    

    return bp         

def boundary_neumann(node1,node2,g_neumann,quad_degre):
    """Return 2*1 elemental neumann boundary vector"""
    node1 = node1[0:2]
    node2 = node2[0:2]
    Xp,W = quadrature_segment(quad_degre)
    sigma = np.linalg.norm(np.array(node2) -np.array(node1))
    bp = np.matrix('0.0 0.0').T
    for i in range(quad_degre):
        xm = (1-Xp[i])*np.array(node1)+Xp[i]*np.array(node2)
        bp[0] += W[i]*g_neumann(xm[0],xm[1])*(1-Xp[i])

    for i in range(quad_degre):
        xm = (1-Xp[i])*np.array(node1)+Xp[i]*np.array(node2)
        bp[1] += W[i]*g_neumann(xm[0],xm[1])*Xp[i]    

    return bp

def boundary_dirichlet_A(elements,A,bord_in_tag):
    Atemp = A.toarray()
    for k in range(1,len(elements)):
        if elements[k][1][0] == bord_in_tag:
            Atemp[elements[k][2][0],:] = 0
            Atemp[elements[k][2][0],elements[k][2][0]] = 1
            Atemp[elements[k][2][1],:] = 0
            Atemp[elements[k][2][1],elements[k][2][1]] = 1
    return csr_matrix(A)

def boundary_elemental_mass_matrix(node1,node2,cst):
    """Return 2*2 elemental matrix fourier-robin condition"""
    sigma = np.linalg.norm(np.array(node2) -np.array(node1))
    sigma = sigma
    Mref = 1./6*np.matrix('2 1;1 2')
    Me = cst*sigma*Mref
    return Me

def assemble(nodes,elements,f,g_neumann,u_dirichlet,k,signe,quad_degre,cst,bord_out_tag,bord_in_tag,interior):
    """assembly sparse Matrix A and second member B """
    rows = []
    cols = []
    vals = []
    B = np.zeros(len(nodes),dtype=complex) # second member
    for k in range(1,len(elements)+1):
        if elements[k][1][0] == interior:
            Me=elemental_mass_matrix(nodes[elements[k][2][0]],nodes[elements[k][2][1]],nodes[elements[k][2][2]],k)
            Dep=elemental_stiffness_matrix(nodes[elements[k][2][0]],nodes[elements[k][2][1]],nodes[elements[k][2][2]],signe)
            if f == 0:
                bp = [0.0,0.0,0.0]
            else:
                bp=elemental_vector(nodes[elements[k][2][0]],nodes[elements[k][2][1]],nodes[elements[k][2][2]],f)    
            for row in range(3):
                glob_row = elements[k][2][row]
                for col in range(3):
                    glob_col = elements[k][2][col]
                    rows.append(glob_row-1)
                    cols.append(glob_col-1)
                    vals.append(Me[row,col]+Dep[row,col])
                B[glob_row-1] += bp[row]         
        elif elements[k][1][0] == bord_in_tag:
            # neumann condition
            x1 = nodes[elements[k][2][0]]
            x2 = nodes[elements[k][2][1]]
            B[elements[k][2][0]] = u_dirichlet(x1[0],x1[1])
            B[elements[k][2][1]] = u_dirichlet(x2[0],x2[1])
        elif elements[k][1][0] == bord_out_tag:
            Me_bord = boundary_elemental_mass_matrix(nodes[elements[k][2][0]],nodes[elements[k][2][1]],cst)
            if g_neumann == 0:
                bp_bord = [0.0,0.0]
            else:    
                bp_bord = boundary_neumann(nodes[elements[k][2][0]],nodes[elements[k][2][1]])
            for row in range(2):
                glob_row = elements[k][2][row]
                for col in range(2):
                    glob_col = elements[k][2][col]
                    rows.append(glob_row-1)
                    cols.append(glob_col-1)
                    vals.append(Me_bord[row,col])
                B[glob_row-1] += bp_bord[row]    
    return  coo_matrix((vals, (rows, cols)), shape=(len(nodes), len(nodes))).tocsr(), B

def solve(A,B):
    u_sol = spsolve(A, B)
    return u_sol