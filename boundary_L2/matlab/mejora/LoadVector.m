function [b1,b2]=LoadVector(mesh,fem,vert_elem)
    
    b1 = zeros(3,1);
    b2 = zeros(3,1);
    x1=mesh.nodes(vert_elem(1),1);
    y1=mesh.nodes(vert_elem(1),2);
    x2=mesh.nodes(vert_elem(2),1);
    y2=mesh.nodes(vert_elem(2),2);
    x3=mesh.nodes(vert_elem(3),1);
    y3=mesh.nodes(vert_elem(3),2);

    b1(1)=fem.integraf1(fem,1,x1,y1,x2,y2,x3,y3);
    b1(2)=fem.integraf1(fem,2,x1,y1,x2,y2,x3,y3);
    b1(3)=fem.integraf1(fem,3,x1,y1,x2,y2,x3,y3);
    b2(1)=fem.integraf2(fem,1,x1,y1,x2,y2,x3,y3);
    b2(2)=fem.integraf2(fem,2,x1,y1,x2,y2,x3,y3);
    b2(3)=fem.integraf2(fem,3,x1,y1,x2,y2,x3,y3);

end
