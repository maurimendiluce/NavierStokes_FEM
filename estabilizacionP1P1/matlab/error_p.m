function error=error_p(nodos,elem,ph)

  ntr=size(elem,1); % cantidad de triangulos
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

    p_h01=(ph(v_elem(1))+ph(v_elem(2)))/2;
    p_h12=(ph(v_elem(2))+ph(v_elem(3)))/2;
    p_h20=(ph(v_elem(3))+ph(v_elem(1)))/2;

    p_m01=p_ex([x01,y01])';
    p_m12=p_ex([x12,y12])';
    p_m20=p_ex([x02,y02])';

    int=int+(abs(delta)/6)*((p_m01-p_h01).^2+(p_m12-p_h12).^2+(p_m20-p_h20).^2);
  end
  error=sqrt(int);
end
