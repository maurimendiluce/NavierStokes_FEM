function D=Dirich_square(e,p,xmin,xmax,ymin,ymax)
  D=[];
  for i=1:length(e)
    nodo=e(i);
    v=p(nodo,:);
    if v(1)==xmin
      D(end+1)=nodo;
    end
    if v(1)==xmax
      D(end+1)=nodo;
    end
    if v(2)==ymin
      D(end+1)=nodo;
    end
    if v(2)==ymax
      D(end+1)=nodo;
    end
  end