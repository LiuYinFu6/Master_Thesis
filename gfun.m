function f = gfun(z,t,y,sig_mul,sig_add)
    f = log(normpdf(z,t,sig_mul).*normpdf(y,exp(z),sig_add));
end