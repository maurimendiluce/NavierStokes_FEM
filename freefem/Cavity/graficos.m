addpath('ffmatlib');

clear all;

addpath('ffmatlib');
[p,b,t]=ffreadmesh('mesh_square.msh');
[vhseq]=ffreaddata('espacio_Vh.txt');
[qhseq]=ffreaddata('espacio_Qh.txt');
[u1,u2]=ffreaddata('solucion.txt');
[dp]=ffreaddata('presion.txt');

%%Mesh
%figure();
%ffpdeplot(p,b,t, ...
%          'Mesh','on', ...
%          'Boundary','on', ...
%          'Title','Mesh');


%%%%%%% 3D Surface
figure();
ffpdeplot(p,b,t, ...
         'XYData',dp,'VhSeq',qhseq,'ZStyle','continuous');%, ...)
          %'Title','Solución ph');

%figure();
%ffpdeplot(p,b,t, ...
%          'XYData',u1,'VhSeq',vhseq,'ZStyle','continuous');%, ...)


%figure();
%subplot(1,2,1)
%ffpdeplot(p,b,t, ...
%          'XYData',u1,'VhSeq',vhseq,'Boundary','off','Contour','on', 'XYStyle', 'off', 'CColor', 'flat', ...
%          'ColorMap','jet', ...
%          'Title','Componente u_1');

%subplot(1,2,2)
%ffpdeplot(p,b,t, ...
%          'XYData',u2,'VhSeq',vhseq,'Boundary','off','Contour','on', 'XYStyle', 'off', 'CColor', 'flat', ...
%          'ColorMap','jet', ...
%          'Title','Componente u_2');

%figure();
%ffpdeplot(p,b,t, ...
%          'FlowData',[u1,u2],'VhSeq',vhseq,'FGridParam',[30,20],'Boundary','on', ...
%          'Title','Solución uh');

axis tight;