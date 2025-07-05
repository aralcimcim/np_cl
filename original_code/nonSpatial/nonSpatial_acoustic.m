% This is to validate the IBP model
% 
% First we test the nonspatial part. 
% 
% The experimental designs that can be used here are two
% acoustic experiments, with 8 sounds (corresponding to the shapes in the
% visual design). Both has 8 sounds that form 4 pairs. In each scene
% (columns in X) 4 sounds are played.
% 
% The design is coded in the nonSpatial_acoustic.mat file
%
% In the joint experiment (X_acJoint) all pairs are shown equal times, as
% in Exp.1,2 in Fiser&Aslin 2001.
%
% In the conditional experiment (X_acCond) there are 2 frequent and 2 rare
% pairs, as in Exp.3 in Fiser&Aslin 2001. 

% noisy-or parameters, testable
lda = 0.9999;
eps = 0.001;

% parameters for new chunks, testable
alpha = 0.005;
Kmax = 1;

% beta prior for r, it's reasonable to leave this U(0,1)
bAlpha= 1;
bBeta = 1;

% code parameters, should show that it's sufficient, higher numbers give
% the same results
wburn = 5;
wsample = 10;
burnIn = 0;

% these don't count for the nonSpatial version
V = [];
sigmaU = 1;
sigmaC = 1;
sigmaV = 1;
phi = 4;


stepNo=100;

%%
tic;
[Zpost_Joint Ypost_Joint Rpost_Joint muCTpost_Joint SigmaCTpost_Joint runTimes_Joint Kmax_Joint wburn_Joint wsample_Joint] = ...
    wood_ibp_learning_frontend_NS(X_acJoint, V, lda, eps, sigmaU, sigmaV, phi, sigmaC, alpha, bAlpha, bBeta, Kmax, wburn, wsample, stepNo, burnIn );
toc
%% 
% converges in 100 steps
% lda = 0.9999; eps = 0.001; alpha = 0.005;
[Zpost_Cond Ypost_Cond Rpost_Cond muCTpost_Cond SigmaCTpost_Cond runTimes_Cond Kmax_Cond wburn_Cond wsample_Cond] = ...
    wood_ibp_learning_frontend_NS(X_acCond, V, lda, eps, sigmaU, sigmaV, phi, sigmaC, alpha, bAlpha, bBeta, Kmax, wburn, wsample, stepNo, burnIn );

%% tests
for i=1:8
    pScene_joint(:,i) = pScene(Xtest_acJoint(:,i), Zpost_Joint, Rpost_Joint, lda, eps);
end

for i=1:12
    pScene_cond(:,i) = pScene(Xtest_acCond(:,i), Zpost_Cond, Rpost_Cond, lda, eps);
end

% compare scenes
beta = 1;
kappa = 0.733;


% joint 
for i=1:4
    pChoice_joint(i) = pChoose_betaKappa(pScene_joint(i), pScene_joint(i+4), beta, kappa);
end

display(pChoice_joint);
mean(pChoice_joint)

% conditional
    pChoice_cond(1) = pChoose_betaKappa(pScene_cond(3), pScene_cond(1), beta, kappa);
    pChoice_cond(2) = pChoose_betaKappa(pScene_cond(4), pScene_cond(2), beta, kappa);
    pChoice_cond(3) = pChoose_betaKappa(pScene_cond(4), pScene_cond(1), beta, kappa);
    pChoice_cond(4) = pChoose_betaKappa(pScene_cond(3), pScene_cond(2), beta, kappa);

display(pChoice_cond);
mean(pChoice_cond)

% single (in the conditional experiment)

for i=1:4
    pChoice_single(i) = pChoose_betaKappa(pScene_cond(i+4), pScene_cond(i+8), beta, kappa);
end
for i=5:7
    pChoice_single(i) = pChoose_betaKappa(pScene_cond(i), pScene_cond(i+5), beta, kappa);
end
    pChoice_single(8) = pChoose_betaKappa(pScene_cond(8), pScene_cond(9), beta, kappa);


display(pChoice_single);
mean(pChoice_single)

%%
