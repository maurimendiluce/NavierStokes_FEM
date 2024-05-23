function stokes

    %%
    % creamos la estructura de datos
    data.mu = 10;
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
        mesh.nodes = p';
        mesh.faces = e';
        mesh.elem = t';
    
        subplot(2,2,1)
        pdemesh(p,e,t)
        title("Malla 2D")
        

        [uh,ph] = solver_stokes(data,mesh,fem);
        subplot(2,2,2)
        pdesurf(mesh.nodes',mesh.elem',ph), shading interp
        title("Solución ph")

        subplot(2,2,3)
        quiver(mesh.nodes(:,1),mesh.nodes(:,2),uh(1:size(mesh.nodes,1)),uh(size(mesh.nodes,1)+1:2*size(mesh.nodes,1)))
        title("Solución uh")

        %startx = -1;
        %starty = -1;
        %subplot(2,2,4)
        %streamline(uh(1:size(mesh.nodes,1)),uh(size(mesh.nodes,1)+1:2*size(mesh.nodes,1)),startx,starty)

        N = N+1;
    end    
   

end