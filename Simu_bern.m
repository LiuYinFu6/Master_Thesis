clc;
clear all;
close all;

load("D:\KULeuven\Master thesis\R\sig_est.mat")
load("D:\KULeuven\Master thesis\R\W_all.mat")
nB = 6;
result = repmat(-999,nB,nB);
save('result.mat','result')
for i = 1:10
    W = W_all(:,i*2-1);
    a = 2;
    b = 2;
    sig1 = sig_est(i,1);
    sig2 = sig_est(i,2);
    [res1,bic] = Estimation_known_lap(W,nB,a,b,sig1,sig2);
    name = ['res_beta34_',num2str(i),'.mat'];
    save(name,'res1');
    name1 = ['bic_beta34_3_',num2str(i),'.mat'];
    save(name1,'bic');
    disp(i)
end


