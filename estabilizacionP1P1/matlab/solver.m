function [wh,ph]=solver(w_inicial,elem,nodos,borde,xmin,xmax,ymin,ymax,medidaOmega,estabiliza,ejemplo,mu)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nn=size(nodos,1);   % numero de nodos

  uh=zeros(2*nn,1);
  ph=zeros(nn,1);

  % Nodos del borde con Cond. de Dirichlet
  NodosD=Dirich_square(borde,nodos,xmin,xmax,ymin,ymax);

  % Nodos libres
  NodosFree=setdiff(1:nn,NodosD);

  [AC,Bx,By,G]=StiffnessMatrix(mu,w_inicial,elem,nodos);

  % Cond. de Dirichlet no homogenea -  u0 dato Dirichlet
  uh(NodosD)=u0(nodos(NodosD,1),nodos(NodosD,2),1,ejemplo); %componente 1
  uh(NodosD+nn)=u0(nodos(NodosD,1),nodos(NodosD,2),2,ejemplo); %componente 2

  if estabiliza==1
    M=[AC zeros(nn) -Bx';zeros(nn) AC -By'; -Bx -By -G];
  else
    M=[AC zeros(nn) -Bx';zeros(nn) AC -By'; -Bx -By -0*G];
  end
  b=LoadVector(elem,nodos,mu,ejemplo);
  
  nodop=NodosFree(1);
  xp=nodos(nodop,1);
  yp=nodos(nodop,2);
  if ejemplo=="Homogeneo"
    ph(nodop)=(2*xp-1)*(2*yp-1);
  end
  if ejemplo=="Cavity"
    ph(nodop)=0;
  end
    
  b=b-M*[uh;ph];  
  K_A1=[AC(NodosFree,NodosFree) zeros(size(NodosFree,2))];
  K_A2=[zeros(size(NodosFree,2)) AC(NodosFree, NodosFree)];
  NPFree=setdiff(1:nn,nodop);
  K_Bx=Bx(NPFree, NodosFree);
  K_By=By(NPFree, NodosFree);
  if estabiliza==1
    K_G=G(NPFree,NPFree);
  else
      K_G=0*G(NPFree,NPFree);
  end
  K=[K_A1 -K_Bx'; K_A2 -K_By'; -K_Bx -K_By -K_G];
  K_bf1=b(NodosFree);
  K_bf2=b(NodosFree+nn);
  K_bf3=b(NPFree+2*nn);
  K_b=[K_bf1;K_bf2;K_bf3];
  wh=K\K_b;

  cant_NodosFree=size(NodosFree,2);
  uh(NodosFree)=wh(1:cant_NodosFree);
  uh(NodosFree+nn)=wh(cant_NodosFree+1:2*cant_NodosFree);
  ph(NPFree)=wh(2*cant_NodosFree+1:2*cant_NodosFree+nn-1);

  if estabiliza==0
    integral=integral_p(nodos,elem,ph);
    ph=ph-integral/medidaOmega;
  end
  if estabiliza==1
    integral_estab=integral_pe(nodos,elem,ph);
    ph=ph-integral_estab/medidaOmega;
  end
  
  wh=uh;
end