include("lib_fem.jl")

#armo la particion y defino los elementos y nodos
N=36
mu=1
ejemplo="Cavity"
nodos,elem,borde = rect_mesh([0,1],[0,1], N, N);

#dato incial
w=zeros(2*length(nodos[:,1]));
#uh,ph=solver(nodos,borde,elem,mu,w,ejemplo)
#metodo de punto fijo
iteraciones=8;
i=0;
while i<iteraciones
    uh,ph=solver(nodos,borde,elem,mu,w,ejemplo);
    w=uh;
    i=i+1;
end