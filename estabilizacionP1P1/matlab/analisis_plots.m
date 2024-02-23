function plot=analisis_plots(uh,ph,nodos,h,estabiliza,ejemplo)

    if ejemplo=="Cavity"
        x=nodos(:,1);
        y=nodos(:,2);
        nn=size(nodos,1);
        [p,e,t]=initmesh('squareR','hmax',h);
        figure(1)
        pdemesh(p,e,t)
        figure(2)
        subplot(121)
        quiver(x,y,uh(1:nn),uh(nn+1:2*nn))
        title("Solucion u_h - Cavity Flow")
        
        subplot(122)
        pdesurf(nodos',elem',ph), shading interp
        if estabiliza=="Si"
            title('Presión p_h estabilizada - Cavity Flow')
        else
            title('Presión p_h no estabilizada - Cavity Flow')
        end
    end
    plot=0;
end
