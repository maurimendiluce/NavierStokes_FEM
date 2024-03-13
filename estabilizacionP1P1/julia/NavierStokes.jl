using LinearAlgebra
using Plots
plotlyjs()
using DelimitedFiles
include("lib_fem.jl")

#importo los datos de la malla
nodos = readdlm("nodes.csv",',')
elem = readdlm("elements.csv",',')
elem=Matrix{Int64}(elem)
elem = elem.+1

n_nodes = size(nodos)[1]
borde = Int64[]
for i=1:n_nodes
    x1,x2 = nodos[i,:]
    if x1==0 
        push!(borde,i)
    elseif x1==1
        push!(borde,i)
    elseif x2==1
        push!(borde,i)
    elseif x2==0
        push!(borde,i)
    end
end

#dato incial
w=zeros(2*length(nodos[:,1]))
ejemplo="Cavity"
mu=1
it = 6 #iteraciones
i=0
uh_1=solver(nodos,borde,elem,mu,w,ejemplo)
uh_2=solver(nodos,borde,elem,mu,uh_1,ejemplo)
uh_3=solver(nodos,borde,elem,mu,uh_2,ejemplo)
uh_4=solver(nodos,borde,elem,mu,uh_3,ejemplo)
uh_5=solver(nodos,borde,elem,mu,uh_4,ejemplo)

x=nodos[:,1]
y=nodos[:,2]
u1=uh_5[1:n_nodes]
u2=uh_5[n_nodes:end]

#quiver(x,y,quiver=(u1,u2))
surface(x,y,u1)
