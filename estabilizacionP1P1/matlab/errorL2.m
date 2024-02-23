function error=errorL2(nodos,elem,uh)
  ntr=size(elem,1); % cantidad de triangulos
  nn=size(nodos,1);   % numero de nodos
  int=0;

  for i = 1 : ntr

    v_elem = elem(i, :);
    v0 = nodos( v_elem(1), : )';
    v1 = nodos( v_elem(2), : )';
    v2 = nodos( v_elem(3), : )';
    x0=v0(1);
    y0=v0(2);
    x1=v1(1);
    y1=v1(2);
    x2=v2(1);
    y2=v2(2);

    delta=(x1-x0)*(y2-y0)+(y0-y1)*(x2-x0);

    x01=(x0+x1)/2;
    x02=(x0+x2)/2;
    x12=(x0+x2)/2;
    y01=(y0+y1)/2;
    y02=(y0+y2)/2;
    y12=(y0+y2)/2;

    uh1_m0=(uh(v_elem(1))+uh(v_elem(2)))/2;
    uh1_m1=(uh(v_elem(2))+uh(v_elem(3)))/2;
    uh1_m2=(uh(v_elem(1))+uh(v_elem(3)))/2;
    uh2_m0=(uh(nn+v_elem(1))+uh(nn+v_elem(2)))/2;
    uh2_m1=(uh(nn+v_elem(2))+uh(nn+v_elem(3)))/2;
    uh2_m2=(uh(nn+v_elem(1))+uh(nn+v_elem(3)))/2;

    [u1_m0,u2_m0]=u_ex([x01,y01]);
    [u1_m1,u2_m1]=u_ex([x12,y12]);
    [u1_m2,u2_m2]=u_ex([x02,y02]);


    dif1_0=uh1_m0-u1_m0;
    dif1_1=uh1_m1-u1_m1;
    dif1_2=uh1_m2-u1_m2;
    dif2_0=uh2_m0-u2_m0;
    dif2_1=uh2_m1-u2_m1;
    dif2_2=uh2_m2-u2_m2;

    int=int+(abs(delta)/6)*(dif1_0^2+dif1_1^2+dif1_2^2+dif2_0^2+dif2_1^2+dif2_2^2);

  end

  error=sqrt(int);
end