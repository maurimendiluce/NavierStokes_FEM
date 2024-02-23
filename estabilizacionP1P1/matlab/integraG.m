function q=integraG(i,j,u1,v1,u2,v2,u3,v3,x1,y1,x2,y2,x3,y3)
  
  if j==1
      delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
  end
  if j==2
      delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
  end
  if j==3
      delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
  end
  x13=(x1+x3)/2;
  y13=(y1+y3)/2;
  x12=(x1+x2)/2;
  y12=(y1+y2)/2;
  x23=(x2+x3)/2;
  y23=(y2+y3)/2;
  proy=1/3;
  a1=(phi(i,x13,y13,u1,v1,u2,v2,u3,v3)-proy)*(phi(j,x13,y13,x1,y1,x2,y2,x3,y3)-proy);
  a2=(phi(i,x12,y12,u1,v1,u2,v2,u3,v3)-proy)*(phi(j,x12,y12,x1,y1,x2,y2,x3,y3)-proy);
  a3=(phi(i,x23,y23,u1,v1,u2,v2,u3,v3)-proy)*(phi(j,x23,y23,x1,y1,x2,y2,x3,y3)-proy);
  q=(a1+a2+a3)*abs(delta)/6;
end
