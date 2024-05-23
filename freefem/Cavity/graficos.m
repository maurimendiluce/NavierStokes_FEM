addpath('ffmatlib');

%mesh
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('mesh_square.msh');
ffpdeplot(p,b,t,'Mesh','on','Boundary','on');
nodos=p(1:2,:)';
elem=t(1:4,:)';

Vh=ffreaddata('espacio_Vh.txt');

ph=ffreaddata('presion.txt');

%figure;
%ffpdeplot(p,b,t,'VhSeq',Vh,'XYData',ph,'Mesh','on','Boundary','on','XLim',[0 1],'YLim',[0 1]);

figure;     
pdesurf(nodos',elem',ph)%, shading interp

[u1,u2]=ffreaddata('velocidad.txt');


%x=nodos(:,1);
%y=nodos(:,2);

figure;
ffpdeplot(p,b,t,'VhSeq',Vh,'XYData',ph,'Mesh','off','Boundary','on', ...
          'XLim',[0 1],'YLim',[0 1],'Contour','on','CColor','b', ...
          'XYStyle','off', 'CGridParam',[150, 150],'ColorBar','off', ...
         'FlowData',[u1,u2],'FGridParam',[24, 24]);