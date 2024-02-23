function error=errorH1(nodos,elem,uh)

  ntr=size(elem,1); % cantidad de triangulos
  nn=size(nodos,1);   % numero de nodos
  int=0;
  uh1=uh(1:nn);
  uh2=uh(nn+1:2*nn);

  grd_bas_fcts = [ -1 -1 ; 1 0 ; 0 1 ]'; 

  for i = 1 : ntr

    v_elem = elem(i, 1:3);
    v1 = nodos( v_elem(1), : )' ;
    v2 = nodos( v_elem(2), : )' ;
    v3 = nodos( v_elem(3), : )' ;

    % We use the cuadrature formula which uses the function values
    % at the midpoint of each side:
    % \int_T  f  \approx  |T| ( f(m12) + f(m23) + f(m31) ) / 3.
    % This formula is exact for quadratic polynomials

    % midpoints of the sides
    m12 = (v1 + v2) / 2;
    m23 = (v2 + v3) / 2;
    m31 = (v3 + v1) / 2;

    % exact gradient at the midpoints of the edges
    aux=grad_u(m12);
    gradu1_12 = aux(1,:);
    aux=grad_u(m23);
    gradu1_23 = aux(1,:);
    aux=grad_u(m31);
    gradu1_31 = aux(1,:);

    aux=grad_u(m12);
    gradu2_12 = aux(2,:);
    aux=grad_u(m23);
    gradu2_23 = aux(2,:);
    aux=grad_u(m31);
    gradu2_31 = aux(2,:);

    % Computation of the gradient of the fe-function uef
    % which is constant on the element.
    % We compute it using the gradients of the basis functions
    % as in the assembly of the system

    B = [ v2-v1 , v3-v1 ];
    area_el = abs(det(B)) / 2;
    
    graduef_1 = B' \ (grd_bas_fcts * uh1(v_elem));
    graduef_2 = B' \ (grd_bas_fcts * uh2(v_elem));

    % differences at the midpoints of the sides
    dif1_12 = graduef_1 - gradu1_12';
    dif1_23 = graduef_1 - gradu1_23';
    dif1_31 = graduef_1 - gradu1_31';

    dif2_12 = graduef_2 - gradu2_12';
    dif2_23 = graduef_2 - gradu2_23';
    dif2_31 = graduef_2 - gradu2_31';

    integral_1 = (dif1_12'*dif1_12 + dif1_23'*dif1_23 + dif1_31'*dif1_31) / 3 ...
        * area_el  ;
    integral_2 = (dif2_12'*dif2_12 + dif2_23'*dif2_23 + dif2_31'*dif2_31) / 3 ...
        * area_el  ;

    int = int + (integral_1+integral_2);

  end

  error = sqrt(int);
end