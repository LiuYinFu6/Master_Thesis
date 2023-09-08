function plotres(a,b,res)
    logA = log(a+b)-log(b);
    logB = log(b);
    s=logB:0.001:(logA+logB);
    x=b:0.001:(a+b);
    plot(x,pdf('beta',(x-b)/a,2,1)./a)
    hold on;
    nB = length(res);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat)
end