function [result,BIC] = Estimation_known_lap(W, nbTerms, a, b, sig_mul,sig_add)
% W: matrix of size nx1 (where n is the sample size), containing the
% observed covariate
% nbTerms: maximum number of terms to be considered in the Bernstein polynomial
% a,b: lower bound and upper bound    

% bound after transformation 
logA = log(a+b)-log(b);
logB = log(b);

% Objects for the output
result = repmat(-999,nbTerms,nbTerms);
BIC = repmat(-999,nbTerms,1);

for nB = 1:nbTerms % nB = 1+order of Bernstein polynomial

    %constraints for the optimization
    beq = 1;
    A = [];
    b = [];
    Aeq = [ones(1,nB)];
    lb = [zeros(nB,1)];
    ub = [ones(nB,1)];
    mynlcon = [];
  
    % starting values for the coefficients of the Bernstein polynomial
    xtemp=unifrnd(0,1,nB,9); 
    xtemp2 = xtemp./repmat(sum(xtemp,1),nB,1);
    x01 = [xtemp2];

    optionsIP = optimoptions(@fmincon,'MaxIter',300,'Algorithm','interior-point','Display','iter','UseParallel',true);
    %optionsAS = optimoptions(@fmincon,'MaxIter',300,'Algorithm','sqp','Display','iter','UseParallel',true);

    g_lik = @(coefsig) likelihood_lap([coefsig;sig_mul;sig_add],nB-1,W,logA,logB);
    [res2,fval2,exitflag2,output2] = fmincon(g_lik,x01(:,1),A,b,Aeq,beq,lb,ub,mynlcon,optionsIP);

    for j = [3,5,7,9]
        [res22,fval22,exitflag22,output22] = fmincon(g_lik,x01(:,j),A,b,Aeq,beq,lb,ub,mynlcon,optionsIP);
        if(fval22<fval2)
            res2=res22;
            fval2=fval22;
            exitflag2=exitflag22;
        end
    end
    result(:,nB) = [res2;zeros(nbTerms-length(res2),1)];
    BIC(nB,:) = 2*fval2+log(length(W))*(nB);
    end
end