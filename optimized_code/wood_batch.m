% parameter definition
lda = 0.999;
eps = 0.01;

sigmaU = 3;
sigmaV = 3;
sigmaC = 0.1;

phi = 4;
bAlpha = 1;
bBeta = 1;

alpha = 0.1;

% loading the inputs
load('./inputs/inputs_2001_exp1.mat', 'trainingdata')
X=trainingdata.shapes;
V=trainingdata.pos;

[N,T] = size(X);
trainingLength = T;

Kmax = 1;
wburn = 5;
wsample = 10;
stepNo = 50;
burnIn = 0;

% loading the lightspeed-master package
addpath(genpath('./lightspeed-master'));

[Zpost, Ypost, Rpost, muCTpost, SigmaCTpost, runTimes, Kmax, wburn, wsample] = wood_ibp_learning_frontend(X(:,1:trainingLength),V(:,1:trainingLength,:),lda,eps,sigmaU,sigmaV,phi,sigmaC,alpha,bAlpha,bBeta, Kmax, wburn, wsample, stepNo, burnIn);
