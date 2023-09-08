function f = g2d(z,t,y,sig_mul,sig_add)
   %f = (normpdf(z,t,sig_mul).*normpdf(y,exp(z),sig_add))*(((sig_add^2*(z-t)-sig_mul^2*(y-exp(z))*exp(z)/(2*sig_mul^2*sig_add^2)))^2-((sig_add^2+2*sig_mul^2*exp(2*z)-sig_mul^2*y)/sig_mul^2*sig_add^2));
   %f = -(sig_add^2+2.*sig_mul^2.*exp(2.*z)-sig_mul^2.*y.*exp(z))./(sig_mul^2.*sig_add^2);
   f = -(1./sig_mul^2)+(y.*exp(z)-2.*exp(2.*z))./(sig_add^2);
end
