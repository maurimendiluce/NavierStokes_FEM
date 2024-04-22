classdef fem

    properties
    end

    methods (Static)
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
    
        
    end
end