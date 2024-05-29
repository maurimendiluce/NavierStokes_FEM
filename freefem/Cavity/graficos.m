addpath('ffmatlib');

clear all;

addpath('ffmatlib');
[p,b,t]=ffreadmesh('mesh_square.msh');
[vhseq]=ffreaddata('espacio_Vh.txt');
[qhseq]=ffreaddata('espacio_Qh.txt');
[u1,u2]=ffreaddata('solucion.txt');
[dp]=ffreaddata('presion.txt');

%%Mesh
figure();
ffpdeplot(p,b,t, ...
          'Mesh','on', ...
          'Boundary','on', ...
          'Title','Mesh');


%%%%%%% 3D Surface
figure();
ffpdeplot(p,b,t, ...
          'XYData',dp,'VhSeq',qhseq,'ZStyle','continuous', ...
          'Title','Solución ph');


figure();
ffpdeplot(p,b,t, ...
          'FlowData',[u1,u2],'VhSeq',vhseq,'FGridParam',[26,30],'Boundary','on', ...
          'Title','Solución uh');

axis tight;