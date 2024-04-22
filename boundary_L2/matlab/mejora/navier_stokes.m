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
mesh.nodes = p';
mesh.faces = e';
mesh.elem = t';

%refinamiento
N=0;
while N <= data.n_ref
    [p,e,t] = refinemesh("squareR",p,e,t,"longest");
    N = N+1;
end    
mesh.nodes = p';
mesh.faces = e';
mesh.elem = t';

figure(1)
pdemesh(p,e,t)




