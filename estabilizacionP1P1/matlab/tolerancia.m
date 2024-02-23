function tol=tolerancia(elem,nodos,w,uh)

  %esto hace mide el error entre w y uh en norma H1 para usar una tolerancia en el metodo de punto fijo
  ntr = size(elem,1);
  nn=size(nodos,1);
  int = 0;

    grd_bas_fcts = [ -1 -1 ; 1 0 ; 0 1 ]';
    
    for i = 1 : ntr

        v_elem = elem(i, 1:3);
        v1 = nodos( v_elem(1), : )' ;
        v2 = nodos( v_elem(2), : )' ;
        v3 = nodos( v_elem(3), : )' ;

        % Computation of the gradient of the fe-function uef
        % which is constant on the element.
        % We compute it using the gradients of the basis functions
        % as in the assembly of the system

        B = [ v2-v1 , v3-v1 ];
        area_el = abs(det(B)) / 2;

        graduh_1 = B' \ (grd_bas_fcts * uh(v_elem));
        graduh_2 = B' \ (grd_bas_fcts * uh(nn+v_elem));
        gradw_1 = B' \ (grd_bas_fcts * w(v_elem));
        gradw_2 = B' \ (grd_bas_fcts * w(nn+v_elem));

        dif1=graduh_1-gradw_1;
        dif2=graduh_2-gradw_2;

        integral = (dif1'*dif1+dif2'*dif2) / 3 ...
            * area_el  ;

        int = int + (integral);

    end

  tol=sqrt(int);
end