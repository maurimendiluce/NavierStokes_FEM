function f1=fun1(x,y,mu,ejemplo)
  if ejemplo=="Homogeneo" %Ejemplo con dato de borde homogeneo
    laplaciano=2*(12*x.^2-12*x+2).*(2*y.^3-3*y.^2+y)+2*(x.^4-2*x.^3+x.^2).*(12*y-6);
    termino_no_lineal=4*(x.^4-2*x.^3+x.^2).*(4*x.^3-6*x.^2+2*x).*(2*y.^3-3*y.^2+y).^2-4*(y.^4-2*y.^3+y.^2).*(2*x.^3-3*x.^2+x).*(x.^4-2*x.^3+x.^2).*(6*y.^2-6*y+1);
    grad_p=2*(2*y-1);
    f1=-mu*laplaciano+termino_no_lineal+grad_p;
  end
  if ejemplo=="Cavity"
    f1=0;% cavity flow
  end