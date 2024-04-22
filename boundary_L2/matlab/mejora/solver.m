function [uh,ph]=solver(u_old,data,mesh,fem)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nn = size(mesh.nodes,1); % numero de nodos
  nt = size(mesh.elem,1);

  uh=zeros(2*nn,1);
  ph=zeros(nn,1);

  % Nodos del borde con Cond. de Dirichlet
  D=mesh.dirich_square(data,mesh);
  NodesD = unique(D);
  % Nodos libres
  NodesFree=setdiff(1:nn,NodesD);

  A = sparse(nn,nn);
  Bx = sparse(nn,nn);
  By = sparse(nn,nn);
  C = sparse(nn,nn);
  G = spare(nn,nn);
  b1 = zeros(nn,1);
  b2 = zeros(nn,1);
    
  for el=1:nt

      vert_elem = mesh.elem(el,:);

      A = A+assamble_laplacian_operator(mesh,fem,vert_elem);
      
      [Bx_aux,By_aux] = assamble_div_operator(mesh,fem,vert_elem);

      Bx = Bx + Bx_aux;
      By = By + By_aux;

      G = G + assamble_estab_operator(mesh,fem,vert_elem);
      
      [b1_aux,b2_aux] = LoadVector(mesh,fem,vert_elem);

      b1 = b1 + b1_aux;
      b2 = b2 + b2_aux;

      
  end

  A = data.mu*A;

  % Cond. de Dirichlet no homogenea -  u0 dato Dirichlet
  uh(NodesD)=fem.u0(mesh.nodes(NodesD,1),mesh.nodes(NodesD,2),1,data); %componente 1
  uh(NodesD+nn)=fem.u0(mesh.nodes(NodesD,1),mesh.nodes(NodesD,2),2,data); %componente 2
        
  M=[A zeros(nn) -Bx';zeros(nn) A -By'; -Bx -By -G];
  b=[b1;b2;zeros(nn,1)];

  nodep = NodesFree(1);
  ph(nodep) = 0;

  b=b-M*[uh;ph];  
  K_A1=[A(NodesFree,NodesFree) zeros(size(NodesFree,2))];
  K_A2=[zeros(size(NodesFree,2)) A(NodesFree, NodesFree)];
  NPFree=setdiff(1:nn,nodop);
  K_Bx=Bx(NPFree, NodesFree);
  K_By=By(NPFree, NodesFree);
  K_G=G(NPFree,NPFree);

  K=[K_A1 -K_Bx'; K_A2 -K_By'; -K_Bx -K_By -K_G];
  K_bf1=b(NodesFree);
  K_bf2=b(NodesFree+nn);
  K_bf3=b(NPFree+2*nn);
  F=[K_bf1;K_bf2;K_bf3];
  u=K\F;

  cant_NodesFree=size(NodesFree,2);
  uh(NodesFree)=u(1:cant_NodesFree);
  uh(NodesFree+nn)=u(cant_NodesFree+1:2*cant_NodesFree);
  ph(NPFree)=u(2*cant_NodesFree+1:2*cant_NodesFree+nn-1);

  medidaOmega = (data.ymax-data.ymin)*(data.xmax-data.xmin);
  integral=fem.integral_p(mesh,fem,ph);
  ph=ph-integral/medidaOmega;
  
  
end
        

