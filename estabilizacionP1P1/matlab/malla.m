%malla
%https://www.mathworks.com/help/pde/ug/pde.pdemodel.generatemesh.html
%https://www.mathworks.com/help/pde/referencelist.html?type=function&category=geometry-and-mesh&s_tid=CRUX_topnav
%https://www.mathworks.com/help/pde/ug/generate-mesh.html

model = createpde;
geometryFromEdges(model,@squareR);
generateMesh(model,"Hmax",0.1,"GeometricOrder","linear");
pdeplot(model)

[p,e,t]=meshToPet(model.Mesh);

