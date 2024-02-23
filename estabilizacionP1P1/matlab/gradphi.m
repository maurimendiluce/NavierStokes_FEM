function [grx,gry]=gradphi(j,x1,y1,x2,y2,x3,y3)
  if j==1
    delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
    grx=(y3-y2)/delta;
    gry=(x2-x3)/delta;
  end
  if j==2
    delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
    grx=(y1-y3)/delta;
    gry=(x3-x1)/delta;
  end
  if j==3
    delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
    grx=(y1-y2)/delta;
    gry=(x2-x1)/delta;
  end
end
