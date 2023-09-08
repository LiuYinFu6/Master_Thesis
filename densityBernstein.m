function f = densityBernstein(x,nB,coef)
% Density of the Bernstein
% nB = order of Bernstein polynomial
% x: the input variable
% coef: coefficients of Bernstein polynomial
if nB==0
    f = coef*betapdf(x,1,1);
else
    k=0:1:nB;
    parAlpha=repmat(k'+1,1,length(x)); 
    parBeta=repmat(nB-k'+1,1,length(x));
    f = reshape(coef,1,length(coef))*betapdf(repmat(x,nB+1,1),parAlpha,parBeta);
end