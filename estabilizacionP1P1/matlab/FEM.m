%%
tic
disp('----------------------------------------------------------------')
disp('Programa de Elementos Finitos para el problema de Navier-Stokes')
disp('Estabilización P1P1 en [0,1]x[0,1]')
disp('----------------------------------------------------------------')
disp('----------------------------------------------------------------')

epsilon=10^-8; %tolerancia
[uh,ph,nodos,borde,elem]=NavierStokes(h,xmin,xmax,ymin,ymax,mu,epsilon,estabiliza,ejemplo);

%analisis_plots(uh,ph,nodos,h,estabiliza,ejemplo)
%%
%PLOTS
x=nodos(:,1);
y=nodos(:,2);
nn=size(nodos,1);
[p,e,t]=initmesh('squareR','hmax',h);
figure(1)
pdemesh(p,e,t)
if ejemplo=="Cavity"
    figure(2)
    subplot(2,2,1)
    pdesurf(nodos',elem',uh(nn+1:2*nn)), shading interp
    subplot(2,2,2)
    quiver(x,y,uh(1:nn),uh(nn+1:2*nn))
    title("Solucion u_h - Cavity Flow")
    subplot(2,2,3)
    pdesurf(nodos',elem',ph), shading interp
    if estabiliza=="Si"
        title('Presión p_h estabilizada - Cavity Flow')
    else
        title('Presión p_h no estabilizada - Cavity Flow')
    end
end
if ejemplo=="Homogeneo"
    figure(2)
    subplot(121)
    quiver(x,y,uh(1:nn),uh(nn+1:2*nn))
    title("Solucion u_h - Dato Dirichlet Homogeneo")
    subplot(122)
    pdesurf(nodos',elem',ph), shading interp
    if estabiliza=="Si"
        title('Presión p_h estabilizada - Dato Dirichlet Homogeneo')
    else
        title('Presión p_h no estabilizada - Dato Dirichlet Homogeneo')
    end
end
%%
%CALCULO de ERROR y ORDEN
if ejemplo=="Homogeneo"
    l=length(vector_h);
    errorH1_u=zeros(1,l);
    error_u=zeros(1,l);
    error_p_L2=zeros(1,l);
    
    for i=1:l
        [uh,ph,nodos,borde,elem]=NavierStokes(vector_h(i),xmin,xmax,ymin,ymax,mu,epsilon,estabiliza,ejemplo);
        errorH1_u(i)=errorH1(nodos,elem,uh);
        errorL2_u(i)=errorL2(nodos,elem,uh);
        errorL2_p(i)=error_p(nodos,elem,ph);
    end
    
    %%creo tablacon datos
    h_max=["1/10";"1/15";"1/20";"1/25";"1/30";"1/35";"1/40";"1/45"];
    errorH1_uh=errorH1_u';
    errorL2_uh=errorL2_u';
    errorL2_ph=errorL2_p';
    T=table(h_max,errorH1_uh,errorL2_uh,errorL2_ph)
    %writetable(T,'errores.txt','Delimiter','\t');
    figure(3)
    subplot(121)
    plot(log(vector_h),log(errorH1_u),'--o')
    hold on
    plot(log(vector_h),log(errorL2_u),'g--v')
    hold on
    plot(log(vector_h),log(vector_h),'k')
    hold on
    plot(log(vector_h),2*log(vector_h),'k:')
    hold on
    title("Orden de convergencia para u_h")
    xlabel('log(h)')
    ylabel('log(||error||)')
    legend('||u-u_h||_H1','||u-u_h||_L2','h','h^2')
    subplot(122)
    plot(log(vector_h),log(errorL2_p),'r--*')
    hold on
    plot(log(vector_h),log(vector_h),'k')
    hold on
    title("Orden de convergencia para p_h")
    xlabel('log(h)')
    ylabel('log(||error||)')
    legend('||p-p_h||_L2','h')
end
%%
disp('Tiempo de ejecución del programa (sec):')
display(toc())