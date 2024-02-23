function b=LoadVector(elem,nodos,mu,ejemplo)
  ntr=size(elem,1);
  nn=size(nodos,1);

  b_f1=zeros(nn,1);
  b_f2=zeros(nn,1);

  for i=1:ntr
    vert_elem=elem(i,1:3);
    x1=nodos(vert_elem(1),1);
    y1=nodos(vert_elem(1),2);
    x2=nodos(vert_elem(2),1);
    y2=nodos(vert_elem(2),2);
    x3=nodos(vert_elem(3),1);
    y3=nodos(vert_elem(3),2);

    b_f1(vert_elem(1))=b_f1(vert_elem(1))+integraf1(1,x1,y1,x2,y2,x3,y3,mu,ejemplo);
    b_f1(vert_elem(2))=b_f1(vert_elem(2))+integraf1(2,x1,y1,x2,y2,x3,y3,mu,ejemplo);
    b_f1(vert_elem(3))=b_f1(vert_elem(3))+integraf1(3,x1,y1,x2,y2,x3,y3,mu,ejemplo);
    b_f2(vert_elem(1))=b_f2(vert_elem(1))+integraf2(1,x1,y1,x2,y2,x3,y3,mu,ejemplo);
    b_f2(vert_elem(2))=b_f2(vert_elem(2))+integraf2(2,x1,y1,x2,y2,x3,y3,mu,ejemplo);
    b_f2(vert_elem(3))=b_f2(vert_elem(3))+integraf2(3,x1,y1,x2,y2,x3,y3,mu,ejemplo);
  end

  b=[b_f1;b_f2;zeros(nn,1)];
end