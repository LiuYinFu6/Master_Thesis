 function f = gInt(t,y,sig_mul,sig_add)
    fun = @(z) g1d(z,t,y,sig_mul,sig_add);
    opts = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','Display','off');
    z0 = fsolve(fun,t,opts);
    %f = exp(gfun(z0,t,y,sig_mul,sig_add))*sqrt(2*pi/abs(g2d(z0,t,y,sig_mul,sig_add)));
    f = normpdf(z0,t,sig_mul).*normpdf(y,exp(z0),sig_add).*sqrt(2*pi./abs(g2d(z0,t,y,sig_mul,sig_add)));
end