from func_mesh import *
import numpy as np

class data:

    """
    define los datos del problema
    """

    def __init__(self,h,mu,xmin,xmax,ymin,ymax,ejemplo):
        self.h = h
        self.mu = mu
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.ejemplo = ejemplo 

datos = data(0.1,0.1,-1,1,-1,1,"Cavity")
malla = create_mesh(datos)
nodes = xnod_from_msh(malla, dim=2)
elem = LaG_from_msh(malla)
borde = dirichlet(nodes,datos)


#print(borde)
#plot de la malla
#plot_msh(malla, '2D')
#plt.show()
        