function y=g1_eps(x,eps)
    
    y = 1-sigma(x)*exp(-(1+x)/eps)-sigma(-x)*exp(-(1-x)/eps);
    
end