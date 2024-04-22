function assamble_nolineal_operator(w,data,mesh)
    
    nn = size(mesh.nodes,1);
    nt = size(mesh.elem,1);
    C = sparse(nn,nn);

    for el=1:nt
        
        M=zeros(3);
        vert_elem=mesh.elem(el,:);
        w1= w(vert_elem);
        w2 = w(nn + vert_elem);
        


        
    end

end

    