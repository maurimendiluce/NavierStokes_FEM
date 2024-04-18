classdef fem

    properties
    end

    methods (Static)
        function D=dirich_square(data,mesh)
            D=[];
            for i=1:length(mesh.edges)
                node=mesh.edges(i);
                v=mesh.points(node,:);
                if v(1)==data.xmin
                    D(end+1)=node;
                end
                if v(1)==data.xmax
                    D(end+1)=node;
                end
                if v(2)==data.ymin
                    D(end+1)=node;
                end
                if v(2)==data.ymax
                    D(end+1)=node;
                end
             end
        end
        
        function f1 = f1(x,y) %fuente

            f1 = 0;
        end

        function f2 = f2(x,y) %fuente

            f2 = 0;
        end
        
        function z=phi(j,x,y,x1,y1,x2,y2,x3,y3) %base de lagrange lineal
    
            if j==1
                delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
                z=(x2-x3)*(y-y3)+(y3-y2)*(x-x3);
                z=z/delta;
            end
            if j==2
                delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
                z=(x3-x1)*(y-y1)+(y1-y3)*(x-x1);
                z=z/delta;
            end
            if j==3
                delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
                z=(x2-x1)*(y-y1)+(y1-y2)*(x-x1);
                z=z/delta;
            end
        end 

        function [grx,gry]=gradphi(j,x1,y1,x2,y2,x3,y3) %gradiente de phi
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

        function q=integra_phi(fem,j,x1,y1,x2,y2,x3,y3)
          if j==1
              delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
          end
          if j==2
              delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
          end
          if j==3
              delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
          end
          a1=fem.phi(j,(x1+x3)/2,(y1+y3)/2,x1,y1,x2,y2,x3,y3);
          a2=fem.phi(j,(x1+x2)/2,(y1+y2)/2,x1,y1,x2,y2,x3,y3);
          a3=fem.phi(j,(x2+x3)/2,(y2+y3)/2,x1,y1,x2,y2,x3,y3);
          q=(a1+a2+a3)*abs(delta)/3;
        end

        function q=integraBx(fem,i,j,u1,v1,u2,v2,u3,v3,x1,y1,x2,y2,x3,y3) %integral Bx
          if j==1
              delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
          end
          if j==2
              delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
          end
          if j==3
              delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
          end
          x13=(x1+x3)/2;
          y13=(y1+y3)/2;
          x12=(x1+x2)/2;
          y12=(y1+y2)/2;
          x23=(x2+x3)/2;
          y23=(y2+y3)/2;
          [grx,gry]=fem.gradphi(j,x1,y1,x2,y2,x3,y3);
          a1=fem.phi(i,x13,y13,u1,v1,u2,v2,u3,v3)*grx;
          a2=fem.phi(i,x12,y12,u1,v1,u2,v2,u3,v3)*grx;
          a3=fem.phi(i,x23,y23,u1,v1,u2,v2,u3,v3)*grx;
          q=(a1+a2+a3)*abs(delta)/6;
        end
        
        function q=integraBy(fem,i,j,u1,v1,u2,v2,u3,v3,x1,y1,x2,y2,x3,y3) %integral By
          if j==1
              delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
          end
          if j==2
              delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
          end
          if j==3
              delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
          end
          x13=(x1+x3)/2;
          y13=(y1+y3)/2;
          x12=(x1+x2)/2;
          y12=(y1+y2)/2;
          x23=(x2+x3)/2;
          y23=(y2+y3)/2;
          [grx,gry]=fem.gradphi(j,x1,y1,x2,y2,x3,y3);
          a1=fem.phi(i,x13,y13,u1,v1,u2,v2,u3,v3)*gry;
          a2=fem.phi(i,x12,y12,u1,v1,u2,v2,u3,v3)*gry;
          a3=fem.phi(i,x23,y23,u1,v1,u2,v2,u3,v3)*gry;
          q=(a1+a2+a3)*abs(delta)/6;
        end
        
        function q=integraf1(fem,j,x1,y1,x2,y2,x3,y3) %integral f1
          if j==1
              delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
          end
          if j==2
              delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
          end
          if j==3
              delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
          end
          a1=fem.f1((x1+x3)/2,(y1+y3)/2)*fem.phi(j,(x1+x3)/2,(y1+y3)/2,x1,y1,x2,y2,x3,y3);
          a2=fem.f1((x1+x2)/2,(y1+y2)/2)*fem.phi(j,(x1+x2)/2,(y1+y2)/2,x1,y1,x2,y2,x3,y3);
          a3=fem.f1((x2+x3)/2,(y2+y3)/2)*fem.phi(j,(x2+x3)/2,(y2+y3)/2,x1,y1,x2,y2,x3,y3);
          q=(a1+a2+a3)*abs(delta)/6;
        end
 
        function q=integraf2(fem,j,x1,y1,x2,y2,x3,y3) %integral f2
          if j==1
              delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
          end
          if j==2
              delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
          end
          if j==3
              delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
          end
          a1=fem.f2((x1+x3)/2,(y1+y3)/2)*fem.phi(j,(x1+x3)/2,(y1+y3)/2,x1,y1,x2,y2,x3,y3);
          a2=fem.f2((x1+x2)/2,(y1+y2)/2)*fem.phi(j,(x1+x2)/2,(y1+y2)/2,x1,y1,x2,y2,x3,y3);
          a3=fem.f2((x2+x3)/2,(y2+y3)/2)*fem.phi(j,(x2+x3)/2,(y2+y3)/2,x1,y1,x2,y2,x3,y3);
          q=(a1+a2+a3)*abs(delta)/6;
        end

        function q=integraG(fem,i,j,u1,v1,u2,v2,u3,v3,x1,y1,x2,y2,x3,y3)
          if j==1
              delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
          end
          if j==2
              delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
          end
          if j==3
              delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
          end
          x13=(x1+x3)/2;
          y13=(y1+y3)/2;
          x12=(x1+x2)/2;
          y12=(y1+y2)/2;
          x23=(x2+x3)/2;
          y23=(y2+y3)/2;
          proy=1/3;
          a1=(fem.phi(i,x13,y13,u1,v1,u2,v2,u3,v3)-proy)*(fem.phi(j,x13,y13,x1,y1,x2,y2,x3,y3)-proy);
          a2=(fem.phi(i,x12,y12,u1,v1,u2,v2,u3,v3)-proy)*(fem.phi(j,x12,y12,x1,y1,x2,y2,x3,y3)-proy);
          a3=(fem.phi(i,x23,y23,u1,v1,u2,v2,u3,v3)-proy)*(fem.phi(j,x23,y23,x1,y1,x2,y2,x3,y3)-proy);
          q=(a1+a2+a3)*abs(delta)/6;
        end

        function int=integral_p(mesh,fem,ph)
          int=0;
          ntr=size(mesh.triang,1);
          for i=1:ntr
            vert_elem=mesh.triang(i,1:3);
            x1=mesh.points(vert_elem(1),1);
            y1=mesh.points(vert_elem(1),2);
            x2=mesh.points(vert_elem(2),1);
            y2=mesh.points(vert_elem(2),2);
            x3=mesh.points(vert_elem(3),1);
            y3=mesh.points(vert_elem(3),2);
        
            int= int + ph(vert_elem(1))*fem.integra_phi(fem,1,x1,y1,x2,y2,x3,y3)+ph(vert_elem(2))*fem.integra_phi(fem,2,x1,y1,x2,y2,x3,y3)+...
                      ph(vert_elem(3))*fem.integra_phi(fem,3,x1,y1,x2,y2,x3,y3);
          end
        end

        function u=u0(x,y,coord,data)
          %coord indica la componente de u0
          l=length(x);
       
          if data.ejemplo=="Cavity" %dato de borde para cavity
            u=zeros(1,l);
            if coord==1
              for j=1:l
                if y(j)==1
                  u(j)=1; 
                end
                if y(j)==-1
                  u(j)=0; 
                end
                if x(j)==1
                  u(j)=0; 
                end
                if x(j)==-1
                  u(j)=0; 
                end
              end
            end
          end 
        end
       
        function tol=tolerancia(mesh,w,uh)

          %esto hace mide el error entre w y uh en norma H1 para usar una tolerancia en el metodo de punto fijo
          ntr = size(mesh.triang,1);
          nn=size(mesh.points,1);
          int = 0;
        
            grd_bas_fcts = [ -1 -1 ; 1 0 ; 0 1 ]';
            
            for i = 1 : ntr
        
                v_elem = mesh.triang(i, 1:3);
                v1 = mesh.points( v_elem(1), : )' ;
                v2 = mesh.points( v_elem(2), : )' ;
                v3 = mesh.points( v_elem(3), : )' ;
        
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

        function b=LoadVector(mesh,fem)
          ntr=size(mesh.triang,1);
          nn=size(mesh.points,1);
        
          b_f1=zeros(nn,1);
          b_f2=zeros(nn,1);
        
          for i=1:ntr
            vert_elem=mesh.triang(i,1:3);
            x1=mesh.points(vert_elem(1),1);
            y1=mesh.points(vert_elem(1),2);
            x2=mesh.points(vert_elem(2),1);
            y2=mesh.points(vert_elem(2),2);
            x3=mesh.points(vert_elem(3),1);
            y3=mesh.points(vert_elem(3),2);
        
            b_f1(vert_elem(1))=b_f1(vert_elem(1))+fem.integraf1(fem,1,x1,y1,x2,y2,x3,y3);
            b_f1(vert_elem(2))=b_f1(vert_elem(2))+fem.integraf1(fem,2,x1,y1,x2,y2,x3,y3);
            b_f1(vert_elem(3))=b_f1(vert_elem(3))+fem.integraf1(fem,3,x1,y1,x2,y2,x3,y3);
            b_f2(vert_elem(1))=b_f2(vert_elem(1))+fem.integraf2(fem,1,x1,y1,x2,y2,x3,y3);
            b_f2(vert_elem(2))=b_f2(vert_elem(2))+fem.integraf2(fem,2,x1,y1,x2,y2,x3,y3);
            b_f2(vert_elem(3))=b_f2(vert_elem(3))+fem.integraf2(fem,3,x1,y1,x2,y2,x3,y3);
          end
        
          b=[b_f1;b_f2;zeros(nn,1)];
        end

        function [AC,Bx,By,G]=StiffnessMatrix(data,mesh,fem,w)

          ntr=size(mesh.triang,1); % cantidad de triangulos
          nn=size(mesh.points,1);   % numero de nodos
          A=zeros(nn);
          C=zeros(nn);
          Bx=zeros(nn);
          By=zeros(nn);
          G=zeros(nn);
        
          for t=1:ntr
        
            vert_elem=mesh.triang(t,1:3);
            x1=mesh.points(vert_elem(1),1);
            y1=mesh.points(vert_elem(1),2);
            x2=mesh.points(vert_elem(2),1);
            y2=mesh.points(vert_elem(2),2);
            x3=mesh.points(vert_elem(3),1);
            y3=mesh.points(vert_elem(3),2);
        
            %OPERADOR A
        
            % calculo de los gradientes de beta
            auxg=zeros(2,3);
            [grx,gry]=fem.gradphi(1,x1,y1,x2,y2,x3,y3); %grad beta_1
            auxg(:,1)=[grx;gry];
            [grx,gry]=fem.gradphi(2,x1,y1,x2,y2,x3,y3); %grad beta_2
            auxg(:,2)=[grx;gry];
            [grx,gry]=fem.gradphi(3,x1,y1,x2,y2,x3,y3); %grad beta_3
            auxg(:,3)=[grx;gry];
            delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
        
            %calculo de la matriz A
            aux_a=(abs(delta)/2)*(auxg'*auxg);
            A(vert_elem,vert_elem)=A(vert_elem,vert_elem)+aux_a;
        
        
            %OPERADOR Bx
        
            auxBx=zeros(3);
        
            auxBx(1,1)=fem.integraBx(fem,1,1,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta1 por beta1_x
            auxBx(2,1)=fem.integraBx(fem,2,1,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta2 por beta1_x
            auxBx(3,1)=fem.integraBx(fem,3,1,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta3 por beta1_x
            auxBx(1,2)=fem.integraBx(fem,1,2,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); 
            auxBx(2,2)=fem.integraBx(fem,2,2,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); 
            auxBx(3,2)=fem.integraBx(fem,3,2,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); 
            auxBx(1,3)=fem.integraBx(fem,1,3,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
            auxBx(2,3)=fem.integraBx(fem,2,3,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
            auxBx(3,3)=fem.integraBx(fem,3,3,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
        
            Bx(vert_elem,vert_elem)=Bx(vert_elem,vert_elem)+auxBx;
        
        
            %OPERADOR By
        
            auxBy=zeros(3);
            for j=1:3
                for i=1:3
                    auxBy(i,j)=fem.integraBy(fem,i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
                end
            end
        
        
            By(vert_elem,vert_elem)=By(vert_elem,vert_elem)+auxBy;
        
            %OPERADOR C
        
            w1_elem=w(vert_elem);
            w2_elem=w(nn+vert_elem);
            w1_12=(1/2)*(w1_elem(1)+w1_elem(2));
            w1_13=(1/2)*(w1_elem(1)+w1_elem(3));
            w1_23=(1/2)*(w1_elem(2)+w1_elem(3));
            w2_12=(1/2)*(w2_elem(1)+w2_elem(2));
            w2_13=(1/2)*(w2_elem(1)+w2_elem(3));
            w2_23=(1/2)*(w2_elem(2)+w2_elem(3));
            w_12=[w1_12,w2_12];
            w_13=[w1_13,w2_13];
            w_23=[w1_23,w2_23];
        
        
            % calculo de los gradientes de beta
            [grad_beta1_x,grad_beta1_y]=fem.gradphi(1,x1,y1,x2,y2,x3,y3); %grad beta_1
            grad_beta1= [grad_beta1_x,grad_beta1_y]';
            [grad_beta2_x,grad_beta2_y]=fem.gradphi(2,x1,y1,x2,y2,x3,y3); %grad beta_2
            grad_beta2=[grad_beta2_x,grad_beta2_y]';
            [grad_beta3_x,grad_beta3_y]=fem.gradphi(3,x1,y1,x2,y2,x3,y3); %grad beta_3
            grad_beta3=[grad_beta3_x,grad_beta3_y]';
           
            delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
        
            aux=zeros(3);
        
            aux(1,1)=(abs(delta)/6)*(w_12+w_13)*grad_beta1;
            aux(1,2)=(abs(delta)/6)*(w_12+w_13)*grad_beta2;
            aux(1,3)=(abs(delta)/6)*(w_12+w_13)*grad_beta3;
            aux(2,1)=(abs(delta)/6)*(w_12+w_23)*grad_beta1;
            aux(2,2)=(abs(delta)/6)*(w_12+w_23)*grad_beta2;
            aux(2,3)=(abs(delta)/6)*(w_12+w_23)*grad_beta3;
            aux(3,1)=(abs(delta)/6)*(w_13+w_23)*grad_beta1;
            aux(3,2)=(abs(delta)/6)*(w_13+w_23)*grad_beta2;
            aux(3,3)=(abs(delta)/6)*(w_13+w_23)*grad_beta3;
        
            C(vert_elem,vert_elem)=C(vert_elem,vert_elem)+aux;
        
            %OPERADOR G
            auxG=zeros(3);
            for i=1:3
                for j=1:3
                    auxG(i,j)=fem.integraG(fem,i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
                end
            end
        
            G(vert_elem,vert_elem)=G(vert_elem,vert_elem)+auxG;
        
          end
        
          A=data.mu*A;
          AC=A+C;
        end

    end
end