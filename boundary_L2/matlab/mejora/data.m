classdef data
    %Datos del problema
    properties
        mu  %cte de viscosidad
        xmin
        xmax
        ymin
        ymax
        domain %geometria del dominio
        example %ejemplo: cavity orig o cavity reg
        n_ref %cantidad de veces que queremos refinar
    end

    methods (Static)
                
        function f1 = f1(x,y) %fuente

            f1 = 0;
        end

        function f2 = f2(x,y) %fuente

            f2 = 0;
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


    end
end