import numpy as np
import matplotlib.pyplot as plt
from fenics import *

n = 16
mesh = RectangleMesh(Point(-1,-1), Point(1, 1), n, n)

V= FunctionSpace(mesh,'P',1)

u_D =Expression('1+x[0]*x[0] +2*x[1]*x[1]',degree=2)

def boundary(x, on_boundary): 
    return on_boundary                 

bc= DirichletBC(V,u_D,boundary)

u= TrialFunction(V) 
v= TestFunction(V) 
f= Constant(-6.0) 
a= dot(grad(u),grad(v))*dx 
L= f*v*dx


u= Function(V)

solve(a ==L, u, bc)

plot(u)
#plot(mesh)
plt.show()
