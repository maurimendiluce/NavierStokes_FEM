classdef fun_data

    properties
    end

    methods (Static)
        function D=dirich_square(data,mesh)
            D=[];
            for i=1:length(mesh.faces)
                node=mesh.faces(i);
                v=mesh.nodes(node,:);
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
            D=unique(D);
        end

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