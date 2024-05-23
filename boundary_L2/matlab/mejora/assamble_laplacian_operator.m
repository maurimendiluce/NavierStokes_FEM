function A=assamble_laplacian_operator(mesh,fem,vert_elem)

 
    x1=mesh.nodes(vert_elem(1),1);
    y1=mesh.nodes(vert_elem(1),2);
    x2=mesh.nodes(vert_elem(2),1);
    y2=mesh.nodes(vert_elem(2),2);
    x3=mesh.nodes(vert_elem(3),1);
    y3=mesh.nodes(vert_elem(3),2);

    grad=zeros(2,3);
    [grx,gry]=fem.gradphi(1,x1,y1,x2,y2,x3,y3); %grad beta_1
    grad(:,1)=[grx;gry];
    [grx,gry]=fem.gradphi(2,x1,y1,x2,y2,x3,y3); %grad beta_2
    grad(:,2)=[grx;gry];
    [grx,gry]=fem.gradphi(3,x1,y1,x2,y2,x3,y3); %grad beta_3
    grad(:,3)=[grx;gry];
    delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
        
    %calculo de la matriz A
    A=(abs(delta)/2)*(grad'*grad);
    

end
        