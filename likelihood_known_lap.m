function f = likelihood_known_lap(coefsig,nB,Y,A,B)
% calculate the likelihood function
% coefsig: initial parameters 
% nB: the order of Bernstein
% Y: observations
% A,B bound
    fy = @(y,theta,sig_mul,sig_add,nB,A,B) integral(@(t) gInt(t,y,sig_mul,sig_add).*densityBernstein((t-B)/A,nB,theta)./A,B,A+B,'ArrayValued',true);
    theta = coefsig(1:(nB+1));
    sig_mul = coefsig(nB+2);
    sig_add = coefsig(nB+3);
    fy_res=fy(Y,theta,sig_mul,sig_add,nB,A,B);
    f = -sum(log(fy_res));
    if f == Inf
        f = 99999; 
    end
    %f=fy_res;
end