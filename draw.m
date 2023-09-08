load('result_beta12.mat')
load('bic_beta12.mat')
[am,index]=min(bic_res,[],1);
coef = repmat(-999,50,6);
a=2;b=2;
logA = log(a+b)-log(b);
logB = log(b);
s=logB:0.001:(logA+logB);
x=b:0.001:(a+b);
subplot(3,2,1)
for i=1:50
    coef(i,:) = result(:,index(i),i);
    nB = index(i);
    res = coef(i,:);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat,'Color',[0.3010 0.7450 0.9330],'LineWidth',1.2)
    hold on;
end
plot(x,pdf('beta',(x-b)/a,1,2)./a,'k','LineWidth',2)
title('2*Beta(1,2)+2','FontSize',15)
ylabel('Density') 

subplot(3,2,2)
load('result_beta75.mat')
load('bic_beta75.mat')
result_sum = result;
bic_res_sum = bic_res;
load('result_beta75_2.mat')
load('bic_beta75_2.mat')
result_sum(:,:,12:50) = result(:,:,12:50);
bic_res_sum(:,12:50) = bic_res(:,12:50);
result = result_sum;
bic_res = bic_res_sum;
[am,index]=min(bic_res,[],1);
coef = repmat(-999,50,6);
a=2;b=2;
logA = log(a+b)-log(b);
logB = log(b);
s=logB:0.001:(logA+logB);
x=b:0.001:(a+b);
hold on;
for i=1:50
    coef(i,:) = result(:,index(i),i);
    nB = index(i);
    res = coef(i,:);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat,'Color',[0.3010 0.7450 0.9330],'LineWidth',1.2)
    hold on;
end
plot(x,pdf('beta',(x-b)/a,0.7,0.5)./a,'k','LineWidth',2)
title('2*Beta(0.7,0.5)+2','FontSize',15)
ylabel('Density') 

subplot(3,2,3)
load('result_beta32.mat')
load('bic_beta32.mat')
[am,index]=min(bic_res,[],1);
coef = repmat(-999,50,6);
a=1;b=2;
logA = log(a+b)-log(b);
logB = log(b);
s=logB:0.001:(logA+logB);
x=b:0.001:(a+b);
hold on;
for i=1:50
    coef(i,:) = result(:,index(i),i);
    nB = index(i);
    res = coef(i,:);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat,'Color',[0.3010 0.7450 0.9330],'LineWidth',1.2)
    hold on;
end
plot(x,pdf('beta',(x-b)/a,3,2)./a,'k','LineWidth',2)
title('2*Beta(3,2)+1','FontSize',15)
ylabel('Density') 

load('result_norm3.mat')
load('bic_norm3.mat')
[am,index]=min(bic_res,[],1);
coef = repmat(-999,50,6);
a=2;b=2;
logA = log(a+b)-log(b);
logB = log(b);
s=logB:0.001:(logA+logB);
x=b:0.001:(a+b);
pd = makedist('Normal','mu',3,'sigma',1.5);
t = truncate(pd,2,4);
subplot(3,2,4)
hold on;
for i=1:50
    coef(i,:) = result(:,index(i),i);
    nB = index(i);
    res = coef(i,:);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat,'Color',[0.3010 0.7450 0.9330],'LineWidth',1.2)
    hold on;
end
plot(x,pdf(t,x),'k','LineWidth',2)
title('N(3,1.5,2,4)','FontSize',15)
ylabel('Density') 

subplot(3,2,5)
load('result_norm2.mat')
load('bic_norm2.mat')
[am,index]=min(bic_res,[],1);
coef = repmat(-999,50,6);
a=2;b=3;
logA = log(a+b)-log(b);
logB = log(b);
s=logB:0.001:(logA+logB);
x=b:0.001:(a+b);
pd = makedist('Normal','mu',4,'sigma',0.5);
t = truncate(pd,2,5);
hold on;
for i=1:50
    coef(i,:) = result(:,index(i),i);
    nB = index(i);
    res = coef(i,:);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat,'Color',[0.3010 0.7450 0.9330],'LineWidth',1.2)
    hold on;
end
plot(x,pdf(t,x),'k','LineWidth',2)
title('N(4,0.5,3,5)','FontSize',15)
ylabel('Density') 

subplot(3,2,6)
load('result_norm1_1.mat')
load('bic_norm1_1.mat')
result_sum = result;
bic_res_sum = bic_res;
load('result_norm10.mat')
load('bic_norm10.mat')
result_sum(:,:,12:50) = result(:,:,12:50);
bic_res_sum(:,12:50) = bic_res(:,12:50);
result = result_sum;
bic_res = bic_res_sum;
[am,index]=min(bic_res,[],1);
coef = repmat(-999,50,6);
a=2;b=2;
logA = log(a+b)-log(b);
logB = log(b);
s=logB:0.001:(logA+logB);
x=b:0.001:(a+b);
pd = makedist('Normal','mu',3,'sigma',2);
t = truncate(pd,2,4);
hold on;
for i=1:50
    coef(i,:) = result(:,index(i),i);
    nB = index(i);
    res = coef(i,:);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
    %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
    hold on;
    plot(x,x_hat,'Color',[0.3010 0.7450 0.9330],'LineWidth',1.2)
    hold on;
end
plot(x,pdf(t,x),'k','LineWidth',2)
title('N(3,2,2,4)','FontSize',15)
ylabel('Density') 


% load('result_beta34.mat')
% load('bic_beta34.mat')
% [am,index]=min(bic_res,[],1);
% coef = repmat(-999,50,6);
% a=2;b=2;
% logA = log(a+b)-log(b);
% logB = log(b);
% s=logB:0.001:(logA+logB);
% x=b:0.001:(a+b);
% plot(x,pdf('beta',(x-b)/a,3,4)./a,'k','LineWidth',3)
% hold on;
% for i=1:50
%     coef(i,:) = result(:,index(i),i);
%     nB = index(i);
%     res = coef(i,:);
%     k=0:1:(nB-1);
%     parAlpha=repmat(k'+1,1,length(x));
%     parBeta=repmat(nB-1-k'+1,1,length(x));
%     x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((log(x)-logB)./logA,nB,1),parAlpha,parBeta).*repmat((1./(logA.*x)),nB,1));
%     %s_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((exp(s)-b)./logA,nB,1),parAlpha,parBeta).*repmat(exp(s).1./(a),nB,1));
%     hold on;
%     plot(x,x_hat,'-b')
%     hold on;
% end
