function f = g1d(z,t,y,sig_mul,sig_add)
    f = sig_add^2.*(z-t)-sig_mul^2.*(y-exp(z)).*exp(z);
end

