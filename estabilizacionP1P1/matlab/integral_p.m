function int=integral_p(nodos,elem,ph)
  int=0;
  ntr=size(elem,1);
  for i=1:ntr
    vert_elem=elem(i,1:3);
    x1=nodos(vert_elem(1),1);
    y1=nodos(vert_elem(1),2);
    x2=nodos(vert_elem(2),1);
    y2=nodos(vert_elem(2),2);
    x3=nodos(vert_elem(3),1);
    y3=nodos(vert_elem(3),2);

    int= int + ph(vert_elem(1))*integra_phi(1,x1,y1,x2,y2,x3,y3)+ph(vert_elem(2))*integra_phi(2,x1,y1,x2,y2,x3,y3)+...
              ph(vert_elem(3))*integra_phi(3,x1,y1,x2,y2,x3,y3);
  end
end