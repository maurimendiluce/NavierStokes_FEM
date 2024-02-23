function phi=phi(j,x,y,x1,y1,x2,y2,x3,y3)
    
    if j==1
        delta=(x2-x3)*(y1-y3)+(y3-y2)*(x1-x3);
        phi=(x2-x3)*(y-y3)+(y3-y2)*(x-x3);
        phi=phi/delta;
    end
    if j==2
        delta=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1);
        phi=(x3-x1)*(y-y1)+(y1-y3)*(x-x1);
        phi=phi/delta;
    end
    if j==3
        delta=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1);
        phi=(x2-x1)*(y-y1)+(y1-y2)*(x-x1);
        phi=phi/delta;
    end
end
