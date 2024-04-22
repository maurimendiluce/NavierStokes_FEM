function [Bx,By]=assamble_div_operator(mesh,fem,vert_elem) 
      
    Bx = zeros(3);
    By = zeros(3);
    
    x1=mesh.nodes(vert_elem(1),1);
    y1=mesh.nodes(vert_elem(1),2);
    x2=mesh.nodes(vert_elem(2),1);
    y2=mesh.nodes(vert_elem(2),2);
    x3=mesh.nodes(vert_elem(3),1);
    y3=mesh.nodes(vert_elem(3),2);
    
    for j=1:3
        for i=1:3
            Bx(i,j)=fem.integraBx(fem,i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
            By(i,j)=fem.integraBy(fem,i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
        end
    end

end