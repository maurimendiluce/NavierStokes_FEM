function y=gh(x,h) %esta funcion sirve para el dominio [-1,1]x[-1,1]
    
    if (-1<=x) && (x<=-1+h)
        y=(x+1)/h;
    elseif (1-h<=x) && (x<=1)
        y=(1-x)/h;
    else
        y=1;
    end
end
    