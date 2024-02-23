function f2=fun2(x,y,mu,ejemplo)
  if ejemplo=="Homogeneo" %Ejemplo con dato de borde homogeneo
    laplaciano=-2*(y.^4-2*y.^3+y.^2).*(12*x-6)-2*(12*y.^2-12*y+2).*(2*x.^3-3*x.^2+x);
    termino_no_lineal=-4*(x.^4-2*x.^3+x.^2).*(2*y.^3-3*y.^2+y).*(y.^4-2*y.^3+y.^2).*(6*x.^2-6*x+1)+4*(y.^4-2*y.^3+y.^2).*(2*x.^3-3*x.^2+x).^2.*(4*y.^3-6*y.^2+2*y);
    grad_p=2*(2*x-1);
    f2=-mu*laplaciano+termino_no_lineal+grad_p;
  end
  if ejemplo=="Cavity"
    f2=0;% cavity flow
  end
