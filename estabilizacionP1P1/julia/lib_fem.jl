#using Delaunay
#using LinearAlgebra
#using GLMakie
#using Plots
#plotlyjs()

#######FUNCIONES PARA GENERAR LA GEOMETRIA. EJEMPLO: RECTANGULO##############
#function rect_mesh(xlim,ylim,N,M)
#    delta_x = (xlim[2] - xlim[1])/N;
#    delta_y = (ylim[2] - ylim[1])/M;
#    points=[0 0]
#    for i in range(0,M)
#        for j in range(0,N)
#            points=vcat(points,[xlim[1]+j*delta_x  ylim[1]+i*delta_y]);
#        end
#    end
#    points=points[setdiff(1:end,1),:];
#    mesh = delaunay(points)
#    elements=mesh.simplices
#    nodos=mesh.points
#
#    borde=[0]
#    for i in range(1,size(nodos)[1])
#        x=nodos[i,1]
#        y=nodos[i,2]
#        if x==xlim[1] || x==xlim[2] || y==ylim[1] || y==ylim[2]
#            borde=vcat(borde,i)
#        end
#    end
#    borde=borde[2:end]
#    return nodos,elements,borde
#end

#####FUNCIONES BASES PARA CADA TRIANGULO#########
function phi(j,X,nodos,elem)
    x=X[1]
    y=X[2]
    vertices=nodos[elem,:]
    x1,x2,x3=vertices[:,1]
    y1,y2,y3=vertices[:,2]
    if j==1 
        delta=(x3-x2)*(y1-y2)+(y2-y3)*(x1-x2)
        return ((x3-x2)*(y-y2)+(y2-y3)*(x-x2))/delta
    elseif j==2
        delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1)
        return ((x3-x1)*(y-y1)+(y1-y3)*(x-x1))/delta
    else
        delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
        return ((x2-x1)*(y-y1)+(y1-y2)*(x-x1))/delta 
    end
end

function gradphi(j,nodos,elem)
    vertices=nodos[elem,:]
    x1,x2,x3=vertices[:,1]
    y1,y2,y3=vertices[:,2]
    if j==1
        delta=(x3-x2)*(y1-y2)+(y2-y3)*(x1-x2)
        return [(y2-y3)/delta,(x3-x2)/delta]
    elseif j==2 
        delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1)
        return [(y1-y3)/delta,(x3-x1)/delta]
    else
        delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
        return [(y1-y2)/delta,(x2-x1)/delta]
    end
end

function puntos_med(nodos,elem)
    vertices=nodos[elem,:]
    x1,x2,x3=vertices[:,1]
    y1,y2,y3=vertices[:,2]
    return 0.5*[[x1+x2,y1+y2],[x1+x3,y1+y3],[x2+x3,y2+y3]]
end

function quad(nodos,elem)
    vertices=nodos[elem,:]
    x1,x2,x3=vertices[:,1]
    y1,y2,y3=vertices[:,2]
    return (x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
end

function integra_phi(j,nodos,elem)
    p_med=puntos_med(nodos, elem)
    a1=phi(j, p_med[1], nodos, elem)
    a2=phi(j, p_med[2], nodos, elem)
    a3=phi(j, p_med[3], nodos, elem)
    
    q=abs(quad(nodos, elem)/3)*(a1+a2+a3)
    return q
end

##MATRICES##

function StiffnessA(nodos,elem)
    # Matriz de rigidez del laplaciano.
    
    # Reservamos espacio
    n_nodes = size(nodos)[1]
    n_elem = size(elem)[1]
    A = zeros(n_nodes,n_nodes)
    
        
    # Matriz global, recorriendo los elementos
    for i in range(1,n_elem)
        el=elem[i,:]

        # Pre-calcular matriz local: int_T0 gradphi_i gradphi_j dx
        S = zeros(2,3)
        S[:,1]=gradphi(0, nodos, el)
        S[:,2]=gradphi(1, nodos, el)
        S[:,3]=gradphi(2, nodos, el)
        S=(abs(quad(nodos, el))/2)*(S'*S)
        A[el,el] = A[el,el] + S
    end
    
    return A
end

function integraB(i,j,grad,nodos,elem)
    p_med=puntos_med(nodos, elem)
    if grad==1 #aca miro phi*gradphi_x 
        phi_x=gradphi(j, nodos, elem)[1]
        a1=phi(i, p_med[1], nodos, elem)*phi_x
        a2=phi(i, p_med[2], nodos, elem)*phi_x
        a3=phi(i, p_med[3], nodos, elem)*phi_x
        
        q=(a1+a2+a3)*abs(quad(nodos, elem))/6
        
        return q
    else
        phi_y=gradphi(j, nodos, elem)[1]
        a1=phi(i, p_med[1], nodos, elem)*phi_y
        a2=phi(i, p_med[2], nodos, elem)*phi_y
        a3=phi(i, p_med[3], nodos, elem)*phi_y
        
        q=(a1+a2+a3)*abs(quad(nodos, elem))/6
        
        return q
    end
end

function StiffnessB(nodos,elem)
    # Reservamos espacio
    
    n_nodes = size(nodos)[1]
    n_elem = size(elem)[1]
    Bx = zeros(n_nodes,n_nodes)
    By = zeros(n_nodes,n_nodes)

    for k in range(1,n_elem)
        el=elem[k,:]
        Sx = zeros(3,3)
        Sy = zeros(3,3)
        for i in range(1,3)
            for j in range(1,3)
                Sx[i,j]=integraB(i, j, 1, nodos, el)
                Sy[i,j]=integraB(i, j, 2, nodos, el)
            end
        end
        Bx[el,el] = Bx[el,el] + Sx
        By[el,el] = By[el,el] + Sy
    end
                
    return Bx,By
end

function StiffnessC(nodos,elem,w) #operador no lineal donde w es la u inicial para iterar

    n_nodes = size(nodos)[1]
    n_elem = size(elem)[1]
    C = zeros(n_nodes,n_nodes)
    
    for k in range(1,n_elem)
        el=elem[k,:]
        w1=w[el]
        w2=w[el.+n_nodes]
        
        w_01=0.5*[w1[1]+w1[2],w2[1]+w2[2]]
        w_02=0.5*[w1[1]+w1[3],w2[1]+w2[3]]
        w_12=0.5*[w1[2]+w1[3],w2[2]+w2[3]]
        
        S=zeros(3,3)
        grad = zeros(2,3)
        grad[:,1]=gradphi(1, nodos, el)
        grad[:,2]=gradphi(2, nodos, el)
        grad[:,3]=gradphi(3, nodos, el)
        S[1,:]=abs(quad(nodos, el)/6)*(w_01+w_02)'*grad
        S[2,:]=abs(quad(nodos, el)/6)*(w_01+w_12)'*grad
        S[3,:]=abs(quad(nodos, el)/6)*(w_02+w_12)'*grad
        
        C[el,el] = C[el,el] + S
    end
    return C
end


function integraG(i,j,nodos,elem)
    p_med=puntos_med(nodos, elem)
    proy=1/3
    a1=(phi(i, p_med[1], nodos, elem)-proy)*(phi(j, p_med[1], nodos, elem)-proy)
    a2=(phi(i, p_med[2], nodos, elem)-proy)*(phi(j, p_med[2], nodos, elem)-proy)
    a3=(phi(i, p_med[3], nodos, elem)-proy)*(phi(j, p_med[3], nodos, elem)-proy)
    
    q=(abs(quad(nodos, elem))/6)*(a1+a2+a3)
    
    return q
end

function StiffnessG(nodos,elem)
    n_nodes = size(nodos)[1]
    n_elem = size(elem)[1]
    G = zeros(n_nodes,n_nodes)
    
    for k in range(1,n_elem)
        el=elem[k,:]
        S = zeros(3,3)
        for i in range(1,3)
            for j in range(1,3)
                S[i,j]=integraG(i, j, nodos, el)
            end
        end
        G[el,el] = G[el,el] + S
    end
    return G
end

##LADO DERECHO##

function f(X,mu,ejemplo)
    
    if ejemplo=="Cavity"
        f=[0,0]
    end
    return f
end

function integraf(i,nodos,elem,mu,ejemplo)
    
    p_med=puntos_med(nodos, elem)
    a1=f(p_med[1],mu,ejemplo)*phi(i,p_med[1],nodos,elem)
    a2=f(p_med[2],mu,ejemplo)*phi(i,p_med[2],nodos,elem)
    a3=f(p_med[3],mu,ejemplo)*phi(i,p_med[3],nodos,elem)
    
    q=abs(quad(nodos, elem)/6)*(a1+a2+a3)
    
    return q
end

function LoadVector(nodos,elem,mu,ejemplo)
    n_nodes = size(nodos)[1]
    n_elem = size(elem)[1]
    F = zeros(2*n_nodes)
    
    for k in range(1,n_elem)
        el=elem[k,:]
        for i in range(1,3)
            F[el[i]]=F[el[i]]+integraf(i, nodos, el, mu, ejemplo)[1]
            F[n_nodes+el[i]]=F[n_nodes+el[i]]+integraf(i, nodos, el, mu, ejemplo)[2]
        end
    end
    
    F=vcat(F, zeros(n_nodes))
    return F
end

##SOLVER PARA NAVIER STOKES P1P1

function u0_cavity(X,componente)
    
    l = length(X[:,1])
    u = zeros(l)
    if componente==1
        for j in range(1,l)
            if X[:,2][j]==1
                u[j]=1
            end
            if X[:,2][j]==0
                u[j]=0
            end
            if X[:,1][j]==1
                u[j]=0
            end
            if X[:,1][j]==0
                u[j]=0
            end
        end
        return u
    else
        for j in range(1,l)
            if X[:,2][j]==1
                u[j]=0
            end
            if X[:,2][j]==0
                u[j]=0
            end
            if X[:,1][j]==1
                u[j]=0
            end
            if X[:,1][j]==0
                u[j]=0
            end
        end
        return u
    end
end

function integra_p(nodos,elem,ph)
    n_elem = size(elem)[1]
    integral=0
    for k in range(1,n_elem)
        el=elem[k,:]
        integral=integral+ph[el[1]]*integra_phi(1, nodos, el)+ph[el[2]]*integra_phi(2, nodos, el)+ph[el[3]]*integra_phi(3, nodos, el)
    end
    return integral
end

function solver(nodos,borde,elem,mu,w,ejemplo)

    n_nodes = size(nodos)[1]
    uh = zeros(2*n_nodes)
    ph = zeros(n_nodes)
    
    ##cond de Dirichlet
    X=nodos[borde,:]
    uh[borde]=u0_cavity(X,1) #componente 1
    uh[n_nodes.+borde]=u0_cavity(X,2) #componente 2
    ##
    
    A = mu*StiffnessA(nodos, elem)
    Bx , By = StiffnessB(nodos, elem)
    C = StiffnessC(nodos, elem, w)
    G = StiffnessG(nodos, elem)
    
    ceros = zeros(n_nodes,n_nodes)
    
    M1=[A+C ceros -Bx']
    M2=[ceros A+C -By']
    M3=[-Bx -By -G]
    M=[M1;M2;M3]
    
    F = LoadVector(nodos, elem, mu, ejemplo)

    v = range(1,n_nodes)
    NodosFree = setdiff(v,borde)
    n_NFree = length(NodosFree)
    
    nodo_p = NodosFree[1]
    NpFree = setdiff(v,nodo_p)
    
    
    if ejemplo=="Cavity"
        ph[nodo_p]=0
    end

    F = F-M*[uh;ph]
    AC=A+C
    K_A1 = [AC[NodosFree,NodosFree] zeros(n_NFree,n_NFree)]
    K_A2 = [zeros(n_NFree,n_NFree) AC[NodosFree,NodosFree]]
    K_Bx = Bx[NpFree, NodosFree]
    K_By = By[NpFree, NodosFree]
    K_G = G[NpFree,NpFree]
    K1 = [K_A1 -K_Bx']
    K2 = [K_A2 -K_By']
    K3 = [-K_Bx -K_By -K_G]
    K = [K1;K2;K3]
    
    K_F1 = F[NodosFree]
    K_F2 = F[n_nodes.+NodosFree]
    K_F3 = F[2*n_nodes.+NpFree]
    K_F = [K_F1;K_F2;K_F3]

    wh=K\K_F
    
    uh[NodosFree]=wh[1:n_NFree]
    uh[n_nodes.+NodosFree]=wh[n_NFree+1:2*n_NFree]
    ph[NpFree]=wh[2*n_NFree+1:2*n_NFree+n_nodes-1]
    
    integral_estab=integra_p(nodos,elem,ph)
    ph=ph.-integral_estab #luego habría que dividir por la medida de Omega
    
    return uh,ph
end