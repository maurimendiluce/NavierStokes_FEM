classdef mesh

    properties
        nodes %p
        faces %e
        elem  %t
        dirichlet 
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

    end
end