import numpy as np
from scipy.spatial import Delaunay

#######FUNCIONES PARA GENERAR LA GEOMETRIA. EJEMPLO: RECTANGULO############## 

def rect_mesh(xlim,ylim,I,J):
    x = np.linspace(xlim[0],xlim[1],I)
    y = np.linspace(ylim[0],ylim[1],J)
    X,Y = np.meshgrid(x,y)
    P = np.array([X.flatten(),Y.flatten()]).T
    T = Delaunay(P)
    elem=np.array(T.simplices)
    return T,P,elem

# Función que encuentra el borde de una triangulación
def Boundary(T):
    boundary = np.array([])
    #boundary = set()
    for i in range(len(T.neighbors)):
        for k in range(3):
            if (T.neighbors[i][k] == -1):
                nk1,nk2 = (k+1)%3, (k+2)%3 
                boundary = np.append(boundary,int(T.simplices[i][nk1]))
                boundary = np.append(boundary,int(T.simplices[i][nk2]))
    borde=np.array([])
    for i in range(len(boundary)):
        borde=np.append(borde,int(boundary[i]))
    return borde

#####FUNCIONES BASES PARA CADA TRIANGULO#########

def phi(j,X,nodos,elem):
    x=X[0]
    y=X[1]
    vertices=nodos[elem,:]
    x0,x1,x2=vertices[:,0]
    y0,y1,y2=vertices[:,1]
    if (j==0): 
        delta=(x2-x1)*(y0-y1)+(y1-y2)*(x0-x1)
        return ((x2-x1)*(y-y1)+(y1-y2)*(x-x1))/delta
    elif (j==1): 
        delta=(x2-x0)*(y1-y0)+(y0-y2)*(x1-x0)
        return ((x2-x0)*(y-y0)+(y0-y2)*(x-x0))/delta
    else: 
        delta=(x1-x0)*(y2-y0)+(y0-y1)*(x2-x0)
        return ((x1-x0)*(y-y0)+(y0-y1)*(x-x0))/delta 
    
def gradphi(j,nodos,elem):
    vertices=nodos[elem,:]
    x0,x1,x2=vertices[:,0]
    y0,y1,y2=vertices[:,1]
    if (j==0): 
        delta=(x2-x1)*(y0-y1)+(y1-y2)*(x0-x1)
        return np.array([(y1-y2)/delta,(x2-x1)/delta])
    elif (j==1): 
        delta=(x2-x0)*(y1-y0)+(y0-y2)*(x1-x0)
        return np.array([(y0-y2)/delta,(x2-x0)/delta])
    else: 
        delta=(x1-x0)*(y2-y0)+(y0-y1)*(x2-x0)
        return np.array([(y0-y1)/delta,(x1-x0)/delta])

def puntos_med(nodos,elem):
    vertices=nodos[elem,:]
    x0,x1,x2=vertices[:,0]
    y0,y1,y2=vertices[:,1]
    return 0.5*np.array([[x0+x1,y0+y1],[x0+x2,y0+y2],[x1+x2,y1+y2]])
    
def quad(nodos,elem):
    vertices=nodos[elem,:]
    x0,x1,x2=vertices[:,0]
    y0,y1,y2=vertices[:,1]
    return (x1-x0)*(y2-y0)+(y0-y1)*(x2-x0)


def integra_phi(j,nodos,elem):
    p_med=puntos_med(nodos, elem)
    a0=phi(j, p_med[0], nodos, elem)
    a1=phi(j, p_med[1], nodos, elem)
    a2=phi(j, p_med[2], nodos, elem)
    
    q=abs(quad(nodos, elem)/3)*(a0+a1+a2)
    return q

##MATRICES##

def StiffnessA(T,nodos,elem):
    # Matriz de rigidez del laplaciano.
    
    # Reservamos espacio
    n_nodes = len(T.points)
    n_elem = len(T.simplices)
    A = np.zeros([n_nodes,n_nodes])
    
        
    # Matriz global, recorriendo los elementos
    for i in range(n_elem):
        el=elem[i,:]
        # Índices de los vértices del triángulo i-ésimo (T_i)
        vertex_index = T.simplices[i]
        # Pre-calcular matriz local: int_T0 gradphi_i gradphi_j dx
        S = np.zeros([2,3])
        S[:,0]=gradphi(0, nodos, el)
        S[:,1]=gradphi(1, nodos, el)
        S[:,2]=gradphi(2, nodos, el)
        S=(abs(quad(nodos, el))/2)*np.transpose(S)@S
        A[np.ix_(vertex_index,vertex_index)] = A[np.ix_(vertex_index,vertex_index)] + S
    
    return A


def integraB(i,j,grad,nodos,elem):
    p_med=puntos_med(nodos, elem)
    if grad==0: #aca miro phi*gradphi_x 
        phi_x=gradphi(j, nodos, elem)[0]
        a0=phi(i, p_med[0], nodos, elem)*phi_x
        a1=phi(i, p_med[1], nodos, elem)*phi_x
        a2=phi(i, p_med[2], nodos, elem)*phi_x
        
        q=(a0+a1+a2)*abs(quad(nodos, elem))/6
        
        return q
    
    else:
        phi_y=gradphi(j, nodos, elem)[1]
        a0=phi(i, p_med[0], nodos, elem)*phi_y
        a1=phi(i, p_med[1], nodos, elem)*phi_y
        a2=phi(i, p_med[2], nodos, elem)*phi_y
        
        q=(a0+a1+a2)*abs(quad(nodos, elem))/6
        
        return q
    

def StiffnessB(T,nodos,elem):
    # Reservamos espacio
    
    n_nodes = len(T.points)
    n_elem = len(T.simplices)
    Bx = np.zeros([n_nodes,n_nodes])
    By = np.zeros([n_nodes,n_nodes])

    for k in range(n_elem):
        vertex_index = T.simplices[k]
        el=elem[k,:]
        Sx = np.zeros([3,3])
        Sy = np.zeros([3,3])
        for i in range(3):
            for j in range(3):
                Sx[i,j]=integraB(i, j, 0, nodos, el)
                Sy[i,j]=integraB(i, j, 1, nodos, el)
        Bx[np.ix_(vertex_index,vertex_index)] = Bx[np.ix_(vertex_index,vertex_index)] + Sx
        By[np.ix_(vertex_index,vertex_index)] = By[np.ix_(vertex_index,vertex_index)] + Sy
                
    #B=np.concatenate((Bx,By),axis=0)
    return Bx,By

def StiffnessC(T,nodos,elem,w): #operador no lineal donde w es la u inicial para iterar

    n_nodes = len(T.points)
    n_elem = len(T.simplices)
    C = np.zeros([n_nodes,n_nodes])
    
    for k in range(n_elem):
        vertex_index = T.simplices[k]
        el=elem[k,:]
        w1=w[el]
        w2=w[n_nodes+el]
        
        w_01=0.5*np.array([w1[0]+w1[1],w2[0]+w2[1]])
        w_02=0.5*np.array([w1[0]+w1[2],w2[0]+w2[2]])
        w_12=0.5*np.array([w1[1]+w1[2],w2[1]+w2[2]])
        
        S=np.zeros([3,3])
        grad = np.zeros([2,3])
        grad[:,0]=gradphi(0, nodos, el)
        grad[:,1]=gradphi(1, nodos, el)
        grad[:,2]=gradphi(2, nodos, el)
        S[0,:]=abs(quad(nodos, el)/6)*(w_01+w_02)@grad
        S[1,:]=abs(quad(nodos, el)/6)*(w_01+w_12)@grad
        S[2,:]=abs(quad(nodos, el)/6)*(w_02+w_12)@grad
        
        C[np.ix_(vertex_index,vertex_index)] = C[np.ix_(vertex_index,vertex_index)] + S
        
    return C

def integraG(i,j,nodos,elem):
    p_med=puntos_med(nodos, elem)
    proy=1/3
    a0=(phi(i, p_med[0], nodos, elem)-proy)*(phi(j, p_med[0], nodos, elem)-proy)
    a1=(phi(i, p_med[1], nodos, elem)-proy)*(phi(j, p_med[1], nodos, elem)-proy)
    a2=(phi(i, p_med[2], nodos, elem)-proy)*(phi(j, p_med[2], nodos, elem)-proy)
    
    q=(abs(quad(nodos, elem))/6)*(a0+a1+a2)
    
    return q

def StiffnessG(T,nodos,elem):
    n_nodes = len(T.points)
    n_elem = len(T.simplices)
    G = np.zeros([n_nodes,n_nodes])
    
    for k in range(n_elem):
        vertex_index = T.simplices[k]
        el=elem[k,:]
        S = np.zeros([3,3])
        for i in range(3):
            for j in range(3):
                S[i,j]=integraG(i, j, nodos, el)
        
        G[np.ix_(vertex_index,vertex_index)] = G[np.ix_(vertex_index,vertex_index)] + S
    
    return G


##LADO DERECHO##

def f(X,mu,ejemplo):
    
    if ejemplo=="Cavity":
        f=np.array([0,0])
    
    return f

def integraf(i,nodos,elem,mu,ejemplo):
    
    p_med=puntos_med(nodos, elem)
    a0=f(p_med[0],mu,ejemplo)*phi(i,p_med[0],nodos,elem)
    a1=f(p_med[1],mu,ejemplo)*phi(i,p_med[1],nodos,elem)
    a2=f(p_med[2],mu,ejemplo)*phi(i,p_med[2],nodos,elem)
    
    q=abs(quad(nodos, elem)/6)*(a0+a1+a2)
    
    return q

def LoadVector(T,nodos,elem,mu,ejemplo):
    n_nodes = len(T.points)
    n_elem = len(T.simplices)
    F = np.zeros(2*n_nodes)
    
    for k in range(n_elem):
        el=elem[k,:]
        
        for i in range(3):
            F[el[i]]=F[el[i]]+integraf(i, nodos, el, mu, ejemplo)[0]
            F[n_nodes+el[i]]=F[n_nodes+el[i]]+integraf(i, nodos, el, mu, ejemplo)[1]
    
    F=np.concatenate((F, np.zeros(n_nodes)),axis=0)
    return F

##SOLVER PARA NAVIER STOKES P1P1

def u0_cavity(X,componente):
    
    l = len(X[:,1])
    u = np.zeros(l)
    if componente==0:
        for j in range(l):
            if X[:,1][j]==1:
                u[j]=1
            if X[:,1][j]==0:
                u[j]=0
            if X[:,0][j]==1:
                u[j]=0
            if X[:,0][j]==0:
                u[j]=0
        return u
    else:
        for j in range(l):
            if X[:,1][j]==1:
                u[j]=0
            if X[:,1][j]==0:
                u[j]=0
            if X[:,0][j]==1:
                u[j]=0
            if X[:,0][j]==0:
                u[j]=0
        return u

def integra_p(T,nodos,elem,ph):
    n_elem = len(T.simplices)
    integral=0
    for k in range(n_elem):
        el=elem[k,:]
        integral=integral+ph[el[0]]*integra_phi(0, nodos, el)+ph[el[1]]*integra_phi(1, nodos, el)+ph[el[2]]*integra_phi(2, nodos, el)

    return integral


def solver(T,nodos,borde,elem,mu,w,ejemplo):
    

    n_nodes = len(T.points)
    uh = np.zeros(2*n_nodes)
    ph = np.zeros(n_nodes)
    
    ##cond de Dirichlet
    borde=borde.astype('int')
    X=nodos[borde,:]
    uh[borde]=u0_cavity(X,0) #componente 1
    uh[borde+n_nodes]=u0_cavity(X,1) #componente 2
    ##
    
    A = mu*StiffnessA(T, nodos, elem)
    Bx , By = StiffnessB(T, nodos, elem)
    C = StiffnessC(T, nodos, elem, w)
    G = StiffnessG(T, nodos, elem)
    
    ceros = np.zeros([n_nodes,n_nodes])
    
    M1 = np.concatenate((A+C, ceros, -np.transpose(Bx)),axis=1)
    M2 = np.concatenate((ceros, A+C, -np.transpose(By)),axis=1)
    M3 = np.concatenate((-Bx, -By,-G),axis=1)
    M = np.concatenate((M1,M2,M3),axis=0)
    
    F = LoadVector(T, nodos, elem, mu, ejemplo)
    
    v = np.arange(0,n_nodes,1)
    NodosFree = np.setdiff1d(v,borde)
    n_NFree = len(NodosFree)
    
    nodo_p = NodosFree[1]
    NpFree = np.setdiff1d(v,nodo_p)
    
    
    if ejemplo=="Cavity":
        ph[nodo_p]=0

    F = F-M@np.concatenate((uh,ph),axis=0)
    AC=A+C
    K_A1 = np.concatenate((AC[np.ix_(NodosFree,NodosFree)],np.zeros([n_NFree,n_NFree])),axis=1)
    K_A2 = np.concatenate((np.zeros([n_NFree,n_NFree]),AC[np.ix_(NodosFree,NodosFree)]),axis=1)
    K_Bx = Bx[np.ix_(NpFree, NodosFree)]
    K_By = By[np.ix_(NpFree, NodosFree)]
    K_G = G[np.ix_(NpFree,NpFree)]
    K1 = np.concatenate((K_A1,-np.transpose(K_Bx)),axis=1)
    K2 = np.concatenate((K_A2,-np.transpose(K_By)),axis=1)
    K3 = np.concatenate((-K_Bx,-K_By,-K_G),axis=1)
    
    K = np.concatenate((K1,K2,K3),axis=0)
    
    K_F1 = F[np.ix_(NodosFree)]
    K_F2 = F[np.ix_(NodosFree+n_nodes)]
    K_F3 = F[np.ix_(NpFree+2*n_nodes)]
    K_F = np.concatenate((K_F1,K_F2,K_F3),axis=0)
    
    wh=np.linalg.solve(K,K_F)
    
    uh[NodosFree]=wh[0:n_NFree]
    uh[NodosFree+n_nodes]=wh[n_NFree:2*n_NFree]
    ph[NpFree]=wh[2*n_NFree:2*n_NFree+n_nodes]
    
    integral_estab=integra_p(T,nodos,elem,ph)
    ph=ph-integral_estab #luego habría que dividir por la medida de Omega
    
    return uh,ph
        
    