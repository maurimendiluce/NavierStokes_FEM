function u=u0(x,y,coord,ejemplo)
  %coord indica la componente de u0
  ly=length(y);

  if ejemplo=="Homogeneo" %dato de borde homogeneo u_0=0
    if coord==1
      u=zeros(1,ly);
    end
    if coord==2
      u=zeros(1,ly);
    end
  end

  if ejemplo=="Cavity" %dato de borde para cavity
    u=zeros(1,ly);
    if coord==1
      for j=1:ly
        if y(j)==1
          u(j)=1; %para cavity
        end
        if y(j)==-1
          u(j)=0; % para cavity
        end
        if x(j)==1
          u(j)=0; % para cavity
        end
        if x(j)==-1
          u(j)=0; % para cavity
        end
      end
    else
      for j=1:ly
        if y(j)==1
          u(j)=0; %para cavity
        end
        if y(j)==-1
          u(j)=0;  % para cavity
        end
        if x(j)==1
          u(j)=0; % para cavity
        end
        if x(j)==-1
          u(j)=0; % para cavity
        end
      end
    end
  end

  if ejemplo=="Cavity_reg" 
    u=zeros(1,ly);
    if coord==1
      for j=1:ly
        if y(j)==1
          u(j)=gh(x(j),0.1);%g1_eps(x(j),0.3); %1; %para cavity
        end
        if y(j)==-1
          u(j)=0; % para cavity
        end
        if x(j)==1
          u(j)=0; % para cavity
        end
        if x(j)==-1
          u(j)=0; % para cavity
        end
      end
    else
      for j=1:ly
        if y(j)==1
          u(j)=0; %para cavity
        end
        if y(j)==-1
          u(j)=0;  % para cavity
        end
        if x(j)==1
          u(j)=0; % para cavity
        end
        if x(j)==-1
          u(j)=0; % para cavity
        end
      end
    end
  end
