function q=integraf2(j,x1,y1,x2,y2,x3,y3,mu,ejemplo)
  if j==1
      delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
  end
  if j==2
      delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
  end
  if j==3
      delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
  end
  a1=fun2((x1+x3)/2,(y1+y3)/2,mu,ejemplo)*phi(j,(x1+x3)/2,(y1+y3)/2,x1,y1,x2,y2,x3,y3);
  a2=fun2((x1+x2)/2,(y1+y2)/2,mu,ejemplo)*phi(j,(x1+x2)/2,(y1+y2)/2,x1,y1,x2,y2,x3,y3);
  a3=fun2((x2+x3)/2,(y2+y3)/2,mu,ejemplo)*phi(j,(x2+x3)/2,(y2+y3)/2,x1,y1,x2,y2,x3,y3);
  q=(a1+a2+a3)*abs(delta)/6;
end
