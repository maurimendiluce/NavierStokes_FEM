function G=assamble_estab_operator(mesh,fem,vert_elem)

    x1=mesh.nodes(vert_elem(1),1);
    y1=mesh.nodes(vert_elem(1),2);
    x2=mesh.nodes(vert_elem(2),1);
    y2=mesh.nodes(vert_elem(2),2);
    x3=mesh.nodes(vert_elem(3),1);
    y3=mesh.nodes(vert_elem(3),2);

    G=zeros(3);
    for i=1:3
        for j=1:3
            G(i,j)=fem.integraG(fem,i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
        end
    end
        
end
