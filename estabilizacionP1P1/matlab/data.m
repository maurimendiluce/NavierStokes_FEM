%%  problem  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
h=0.05;
xmin=-1;
xmax=1;
ymin=-1;
ymax=1;
domain="squareR";

mu=1;
estabiliza="Si";
ejemplo="Cavity";
%ejemplo="Homogeneo";

%%
%PARA CALCULO de ERROR y ORDEN
vector_h=[1/10,1/15,1/20,1/25,1/30,1/35,1/40,1/45];
