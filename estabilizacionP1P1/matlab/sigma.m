function y=sigma(x)
    if (-1<=x) && (x<=0.5)
        y = 1;
    elseif (-0.5<=x) && (x<=0.5)
        y = 6*(0.5-x)^5-15*(0.5-x)^4+10*(0.5-x)^3;
    elseif (0.5<=x) && (x<=1)
        y = 0;
    end
end