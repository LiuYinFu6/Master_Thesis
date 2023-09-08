function [result,BIC] = Estimation_unknown_lap(W,nbTerms,a,b)
% W: matrix of size nx1 (where n is the sample size), containing the
% observed covariate
% nbTerms: maximum number of terms to be considered in the Bernstein polynomial
% a,b: lower bound and upper bound 

% bound after transformation 
logA = log(a+b)-log(b);
logB = log(b);

% Objects for the output
result = repmat(-999,nbTerms+2,nbTerms);
BIC = repmat(-999,nbTerms,1);

% starting values for v, the measurement error standard deviation
x01temp = (std(W))*[0.2,0.6,1]; 
[m,n] = meshgrid( x01temp, x01temp' );
[sig0(:,1),sig0(:,2) ] = deal( reshape(m,[],1), reshape(n,[],1) ); 
sig0 = sig0';

for nB = nbTerms % nB = 1+order of Bernstein polynomial

    k=0:1:(nB-1);
    parAlpha = k+1;
    parBeta = nB-1-k+1;
    mean_v = arrayfun(@(x,y) hypergeom(x,y,logA),parAlpha,parAlpha+parBeta);
    exp2_v = arrayfun(@(x,y) hypergeom(x,y,2*logA),parAlpha,parAlpha+parBeta);

    mean_est0 = @(theta) exp(logB).*(reshape(theta,1,length(theta))*mean_v');
    exp2_est0 = @(theta) (exp(2*logB).*reshape(theta,1,length(theta))*(reshape(exp2_v,length(exp2_v),1)));
    
    c = @(x) [mean_est0(x(1:(length(x)-2)))-mean(W);
        mean(W)-mean_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^2/2)-0.1;
        mean_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^2/2)-mean(W)-0.1;
        mean(W)^2-exp2_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^4);
        (mean_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^2/2))^2-exp2_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^4);
        var(W) - exp2_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^4)+(mean_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^2/2))^2+x(length(x))^2-0.1;
        -var(W) + exp2_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^4)-(mean_est0(x(1:(length(x)-2)))*exp(x(length(x)-1)^2/2))^2-x(length(x))^2-0.1
        ];
    ceq = @(x) [];
    nonlcon = @(x) deal(c(x),ceq(x));
    

    %constraints for the optimization
    beq = 1;
%     A = [];
%     b = [];
    A = [mean_v.*exp(logB),0,0];
    b = [mean(W)];
    Aeq = [ones(1,nB),zeros(1,2)];
    %lb = [zeros(nB,1);0.01;0.01];
    lb = [zeros(nB,1);0;0];
    ub = [ones(nB,1);std(W)*2;std(W)];
    mynlcon = [];   

    % starting values for the coefficients of the Bernstein polynomial
    xtemp=unifrnd(0,1,nB,9); 
    xtemp2 = xtemp./repmat(sum(xtemp,1),nB,1);
    x01 = [xtemp2;sig0];

    optionsIP = optimoptions(@fmincon,'MaxIter',300,'Algorithm','interior-point','Display','iter','UseParallel',true);
    %optionsAS = optimoptions(@fmincon,'MaxIter',300,'Algorithm','sqp','Display','iter','UseParallel',true);

    % objective function
    g_lik = @(coefsig) likelihood_lap(coefsig,nB-1,W,logA,logB);
    [res2,fval2,exitflag2,output2] = fmincon(g_lik,x01(:,1),A,b,Aeq,beq,lb,ub,nonlcon,optionsIP);
    %for j = [2:9]
    %for j = [3,5,7,9]
    for j = [3]
        [res22,fval22,exitflag22,output22] = fmincon(g_lik,x01(:,j),A,b,Aeq,beq,lb,ub,nonlcon,optionsIP);
        if(fval22<fval2)
            res2=res22;
            fval2=fval22;
            exitflag2=exitflag22;
        end
    end
    result(:,nB) = [res2;zeros(nbTerms+2-length(res2),1)];
    BIC(nB,:) = 2*fval2+log(length(W))*(nB+2);
    end
end