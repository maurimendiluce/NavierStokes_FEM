classdef navier_stokes

    properties
    end

    methods (Static)
        
        function [wh,ph]=solver(w0,data,mesh,fem)

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          nn=size(mesh.points,1);   % numero de nodos
        
          uh=zeros(2*nn,1);
          ph=zeros(nn,1);
        
          % Nodos del borde con Cond. de Dirichlet
          D=fem.dirich_square(data,mesh);
          NodosD = unique(D);
          % Nodos libres
          NodosFree=setdiff(1:nn,NodosD);
        
          [AC,Bx,By,G]=fem.StiffnessMatrix(data,mesh,fem,w0);
        
          % Cond. de Dirichlet no homogenea -  u0 dato Dirichlet
          uh(NodosD)=fem.u0(mesh.points(NodosD,1),mesh.points(NodosD,2),1,data); %componente 1
          uh(NodosD+nn)=fem.u0(mesh.points(NodosD,1),mesh.points(NodosD,2),2,data); %componente 2
        
          M=[AC zeros(nn) -Bx';zeros(nn) AC -By'; -Bx -By -G];
      
          b=fem.LoadVector(mesh,fem);
          
          nodop=NodosFree(1);
          %xp=mesh.points(nodop,1);
          %yp=mesh.points(nodop,2);
          if data.ejemplo=="Cavity"
            ph(nodop)=0;
          end
            
          b=b-M*[uh;ph];  
          K_A1=[AC(NodosFree,NodosFree) zeros(size(NodosFree,2))];
          K_A2=[zeros(size(NodosFree,2)) AC(NodosFree, NodosFree)];
          NPFree=setdiff(1:nn,nodop);
          K_Bx=Bx(NPFree, NodosFree);
          K_By=By(NPFree, NodosFree);
          K_G=G(NPFree,NPFree);
       
          K=[K_A1 -K_Bx'; K_A2 -K_By'; -K_Bx -K_By -K_G];
          K_bf1=b(NodosFree);
          K_bf2=b(NodosFree+nn);
          K_bf3=b(NPFree+2*nn);
          K_b=[K_bf1;K_bf2;K_bf3];
          wh=K\K_b;
        
          cant_NodosFree=size(NodosFree,2);
          uh(NodosFree)=wh(1:cant_NodosFree);
          uh(NodosFree+nn)=wh(cant_NodosFree+1:2*cant_NodosFree);
          ph(NPFree)=wh(2*cant_NodosFree+1:2*cant_NodosFree+nn-1);
        
          medidaOmega = (data.ymax-data.ymin)*(data.xmax-data.xmin);
          integral=fem.integral_p(mesh,fem,ph);
          ph=ph-integral/medidaOmega;
          
          wh=uh;
        end
        
    end
end