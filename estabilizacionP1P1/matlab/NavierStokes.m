function [uh,ph,nodos,borde,elem]=NavierStokes(h,xmin,xmax,ymin,ymax,mu,epsilon,estabiliza,ejemplo)

  [nodos,e,elem]=initmesh('squareR','hmax',h);
  nodos=nodos';
  e=e';
  borde=e(:,1);
  elem=elem';
  nn=size(nodos,1);
  medidaOmega=(xmax-xmin)*(ymax-ymin);

  w_inicial=zeros(2*nn,1);
  
  tol=10;
  iteraciones=0;
  while tol>epsilon
    [uh,ph]=solver(w_inicial,elem,nodos,borde,xmin,xmax,ymin,ymax,medidaOmega,estabiliza,ejemplo,mu);
    tol=tolerancia(elem,nodos,w_inicial,uh);
    w_inicial=uh;
    iteraciones=iteraciones+1;
  end

  uh=w_inicial;
end