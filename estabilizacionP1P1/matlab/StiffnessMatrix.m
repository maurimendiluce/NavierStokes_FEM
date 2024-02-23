function [AC,Bx,By,G]=StiffnessMatrix(mu,w,elem,nodos)

  ntr=size(elem,1); % cantidad de triangulos
  nn=size(nodos,1);   % numero de nodos
  A=zeros(nn);
  C=zeros(nn);
  Bx=zeros(nn);
  By=zeros(nn);
  G=zeros(nn);

  for t=1:ntr

    vert_elem=elem(t,1:3);
    x1=nodos(vert_elem(1),1);
    y1=nodos(vert_elem(1),2);
    x2=nodos(vert_elem(2),1);
    y2=nodos(vert_elem(2),2);
    x3=nodos(vert_elem(3),1);
    y3=nodos(vert_elem(3),2);

    %OPERADOR A

    % calculo de los gradientes de beta
    auxg=zeros(2,3);
    [grx,gry]=gradphi(1,x1,y1,x2,y2,x3,y3); %grad beta_1
    auxg(:,1)=[grx;gry];
    [grx,gry]=gradphi(2,x1,y1,x2,y2,x3,y3); %grad beta_2
    auxg(:,2)=[grx;gry];
    [grx,gry]=gradphi(3,x1,y1,x2,y2,x3,y3); %grad beta_3
    auxg(:,3)=[grx;gry];
    delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);

    %calculo de la matriz A
    aux_a=(abs(delta)/2)*(auxg'*auxg);
    A(vert_elem,vert_elem)=A(vert_elem,vert_elem)+aux_a;


    %OPERADOR Bx

    auxBx=zeros(3);
    %for j=1:3
    %    for i=1:3
    %        auxBx(i,j)=integraBx(i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
    %    end
    %end

    auxBx(1,1)=integraBx(1,1,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta1 por beta1_x
    auxBx(2,1)=integraBx(2,1,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta2 por beta1_x
    auxBx(3,1)=integraBx(3,1,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta3 por beta1_x
    auxBx(1,2)=integraBx(1,2,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta1 por beta2_x
    auxBx(2,2)=integraBx(2,2,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta2 por beta2_x
    auxBx(3,2)=integraBx(3,2,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3); % integra beta3 por beta2_x
    auxBx(1,3)=integraBx(1,3,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
    auxBx(2,3)=integraBx(2,3,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
    auxBx(3,3)=integraBx(3,3,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);

    Bx(vert_elem,vert_elem)=Bx(vert_elem,vert_elem)+auxBx;


    %OPERADOR By

    auxBy=zeros(3);
    for j=1:3
        for i=1:3
            auxBy(i,j)=integraBy(i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
        end
    end

    %auxBy(1,1)=integraBy(x2,y2,x1,y1,x0,y0,x2,y2,x1,y1,x0,y0); % integra beta1 por beta1_y
    %auxBy(2,1)=integraBy(x0,y0,x2,y2,x1,y1,x2,y2,x1,y1,x0,y0); % integra beta2 por beta1_y
    %auxBy(3,1)=integraBy(x0,y0,x1,y1,x2,y2,x2,y2,x1,y1,x0,y0); % integra beta3 por beta1_y
    %auxBy(1,2)=integraBy(x2,y2,x1,y1,x0,y0,x0,y0,x2,y2,x1,y1); % integra beta1 por beta2_y
    %auxBy(2,2)=integraBy(x0,y0,x2,y2,x1,y1,x0,y0,x2,y2,x1,y1); % integra beta2 por beta2_y
    %auxBy(3,2)=integraBy(x0,y0,x1,y1,x2,y2,x0,y0,x2,y2,x1,y1); % integra beta3 por beta2_y
    %auxBy(1,3)=integraBy(x2,y2,x1,y1,x0,y0,x0,y0,x1,y1,x2,y2);
    %auxBy(2,3)=integraBy(x0,y0,x2,y2,x1,y1,x0,y0,x1,y1,x2,y2);
    %auxBy(3,3)=integraBy(x0,y0,x1,y1,x2,y2,x0,y0,x1,y1,x2,y2);

    By(vert_elem,vert_elem)=By(vert_elem,vert_elem)+auxBy;

    %OPERADOR C

    w1_elem=w(vert_elem);
    w2_elem=w(nn+vert_elem);
    w1_12=(1/2)*(w1_elem(1)+w1_elem(2));
    w1_13=(1/2)*(w1_elem(1)+w1_elem(3));
    w1_23=(1/2)*(w1_elem(2)+w1_elem(3));
    w2_12=(1/2)*(w2_elem(1)+w2_elem(2));
    w2_13=(1/2)*(w2_elem(1)+w2_elem(3));
    w2_23=(1/2)*(w2_elem(2)+w2_elem(3));
    w_12=[w1_12,w2_12];
    w_13=[w1_13,w2_13];
    w_23=[w1_23,w2_23];


    % calculo de los gradientes de beta
    [grad_beta1_x,grad_beta1_y]=gradphi(1,x1,y1,x2,y2,x3,y3); %grad beta_1
    grad_beta1= [grad_beta1_x,grad_beta1_y]';
    [grad_beta2_x,grad_beta2_y]=gradphi(2,x1,y1,x2,y2,x3,y3); %grad beta_2
    grad_beta2=[grad_beta2_x,grad_beta2_y]';
    [grad_beta3_x,grad_beta3_y]=gradphi(3,x1,y1,x2,y2,x3,y3); %grad beta_3
    grad_beta3=[grad_beta3_x,grad_beta3_y]';
    %gradientes=[grad_beta1 grad_beta2 grad_beta3];

    delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);

    aux=zeros(3);

    aux(1,1)=(abs(delta)/6)*(w_12+w_13)*grad_beta1;
    aux(1,2)=(abs(delta)/6)*(w_12+w_13)*grad_beta2;
    aux(1,3)=(abs(delta)/6)*(w_12+w_13)*grad_beta3;
    aux(2,1)=(abs(delta)/6)*(w_12+w_23)*grad_beta1;
    aux(2,2)=(abs(delta)/6)*(w_12+w_23)*grad_beta2;
    aux(2,3)=(abs(delta)/6)*(w_12+w_23)*grad_beta3;
    aux(3,1)=(abs(delta)/6)*(w_13+w_23)*grad_beta1;
    aux(3,2)=(abs(delta)/6)*(w_13+w_23)*grad_beta2;
    aux(3,3)=(abs(delta)/6)*(w_13+w_23)*grad_beta3;

    C(vert_elem,vert_elem)=C(vert_elem,vert_elem)+aux;

    %OPERADOR G
    auxG=zeros(3);
    for i=1:3
        for j=1:3
            auxG(i,j)=integraG(i,j,x1,y1,x2,y2,x3,y3,x1,y1,x2,y2,x3,y3);
        end
    end

    %auxG(1,1)=integraE(x2,y2,x1,y1,x0,y0,x2,y2,x1,y1,x0,y0);
    %auxG(1,2)=integraE(x2,y2,x1,y1,x0,y0,x2,y2,x0,y0,x1,y1);
    %auxG(1,3)=integraE(x0,y0,x1,y1,x2,y2,x2,y2,x1,y1,x0,y0);
    %auxG(2,2)=integraE(x2,y2,x0,y0,x1,y1,x2,y2,x0,y0,x1,y1);
    %auxG(2,3)=integraE(x0,y0,x1,y1,x2,y2,x2,y2,x0,y0,x1,y1);
    %auxG(3,3)=integraE(x0,y0,x1,y1,x2,y2,x0,y0,x1,y1,x2,y2);
    %auxG(2,1)=auxG(1,2);
    %auxG(3,1)=auxG(1,3);
    %auxG(3,2)=auxG(2,3);

    G(vert_elem,vert_elem)=G(vert_elem,vert_elem)+auxG;

  end

  A=mu*A;
  AC=A+C;
end