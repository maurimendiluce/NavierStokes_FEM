%%  problem  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
h=0.03;
xmin=0;
xmax=1;
ymin=0;
ymax=1;
domain="squareR";

mu=1;
estabiliza=1;
ejemplo="Cavity";
%ejemplo="Homogeneo";

%%
%PARA CALCULO de ERROR y ORDEN
vector_h=[1/10,1/15,1/20,1/25,1/30,1/35,1/40,1/45];
