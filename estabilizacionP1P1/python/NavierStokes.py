import matplotlib.pyplot as plt
from lib_fem_mixtos import *


#armo la particion y defino los elementos y nodos
N=24
mu=1
ejemplo="Cavity"
#T,nodos,elem = rect_mesh([0,1],[0,1], N, N)
#borde=Boundary(T)

#print(borde)

#plt.triplot(nodos[:,0], nodos[:,1], T.simplices)
#plt.show()


#dato incial
#w=np.zeros(2*len(T.points))


#metodo de punto fijo
iteraciones=8
i=0
while i<iteraciones:
#    uh,ph=solver(T,nodos,borde,elem,mu,w,ejemplo)
#    w=uh
    i=i+1

#X=nodos[:,0]
#Y=nodos[:,1]
#n=len(nodos[:,0])

#uh_1=uh[0:n]
#uh_2=uh[n:2*n]

#fig1, ax1 = plt.subplots()
#ax1.set_title('Solución u_h')
#Q = ax1.quiver(X, Y, uh_1, uh_2, units='width')
#qk = ax1.quiverkey(Q)

# Creating figure
#fig = plt.figure(figsize =(16, 9)) 
#ax = plt.axes(projection ='3d') 
 
# Creating color map
#my_cmap = plt.get_cmap('hot')
   
# Creating plot
#trisurf = ax.plot_trisurf(X, Y, ph,
#                         cmap = my_cmap)#,
                         #linewidth = 0.2,
                         #antialiased = True,
                         #edgecolor = 'grey') 
#fig.colorbar(trisurf, ax = ax, shrink = 0.5, aspect = 5)
#ax.set_title('Solución $p_h$')
 
# Adding labels
#ax.set_xlabel('X-axis', fontweight ='bold')
#ax.set_ylabel('Y-axis', fontweight ='bold')
#ax.set_zlabel('Z-axis', fontweight ='bold')
     
# show plot
plt.show()