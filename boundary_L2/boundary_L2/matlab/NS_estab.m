%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMA PRINCIPAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% creamos la estructura de datos
data.mu = 1;
data.xmin = -1;
data.xmax = 1;
data.ymin = -1;
data.ymax = 1;
data.domain = 'squareR';
data.ejemplo = 'Cavity';
data.n_ref = input('number of refinements:');


%%
% Creamos la malla inicial
[p,e,t]=initmesh('squareR',"Hmax",inf);
mesh.points = p';
mesh.edges = e';
mesh.triang = t';
mesh.dirichlet = fem.dirich_square(data,mesh);

%refinamiento
N=0;
while N <= data.n_ref
    [p,e,t] = refinemesh("squareR",p,e,t,"longest");
    N = N+1;
end    
mesh.points = p';
mesh.edges = e';
mesh.triang = t';

%%
%Resolvemos usando Newton
nn = size(mesh.points,1);
w0=zeros(2*nn,1);
  
tol=1;
it=0;
while tol>10^-8 %falta poner un criterio de parada en pasos
    [uh,ph] = navier_stokes.solver(w0,data,mesh,fem);
    tol = fem.tolerancia(mesh,w0,uh);
    w0 = uh;
    it = it+1;
end
uh = w0;

%%
%PLOTS
subplot(1,3,1)
pdemesh(p,e,t)

subplot(1,3,2)
pdesurf(mesh.points',mesh.triang',ph), shading interp

subplot(1,3,3)
quiver(mesh.points(:,1),mesh.points(:,2),uh(1:nn),uh(nn+1:2*nn))
