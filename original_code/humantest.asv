function [pTest_SingleEq pTest_Single pTest_Joint pTest_Cond] = humantest(Rtreshold, XTest_Human_SingleEq, VTest_Human_SingleEq, ...
                                                                XTest_Human_Single, VTest_Human_Single, ...
                                                                XTest_Human_Joint, VTest_Human_Joint, ...
                                                                XTest_Human_Cond, VTest_Human_Cond, ...
                                                                Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Vswitch)
if ~exist('Vswitch','var')
    % third parameter does not exist, so default it to something
    Vswitch = 1;
end

                                                            
% test for humans
% ZpostH20 = ZpostH20S;
% RpostH20 = RpostH20S;
% YpostH20 = YpostH20S;

%%
tic
testNum = 4;
pTest_Single = NaN(testNum, 2);
for i = 1:testNum
pTest_Single(i, 1) = pScene_spatial_Cpost(XTest_Human_Single(:, (i-1)*2 + 1), VTest_Human_Single(:,(i-1)*2 + 1,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
pTest_Single(i, 2) = pScene_spatial_Cpost(XTest_Human_Single(:, (i-1)*2 + 2), VTest_Human_Single(:,(i-1)*2 + 2,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch)
end
toc

%%
testNum = 4;
pTest_SingleEq = NaN(testNum, 2);
tic
for i = 1:testNum
    
pTest_SingleEq(i, 1) = pScene_spatial_Cpost(XTest_Human_SingleEq(:, (i-1)*2 + 1), VTest_Human_SingleEq(:,(i-1)*2 + 1,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
pTest_SingleEq(i, 2) = pScene_spatial_Cpost(XTest_Human_SingleEq(:, (i-1)*2 + 2), VTest_Human_SingleEq(:,(i-1)*2 + 2,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
    
end
toc 

%%
testNum = 4;
pTest_Joint = NaN(testNum, 2);
for i = 1:testNum
pTest_Joint(i, 1) = pScene_spatial_Cpost(XTest_Human_Joint(:, (i-1)*2 + 1), VTest_Human_Joint(:,(i-1)*2 + 1,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
pTest_Joint(i, 2) = pScene_spatial_Cpost(XTest_Human_Joint(:, (i-1)*2 + 2), VTest_Human_Joint(:,(i-1)*2 + 2,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
end
%%
testNum = 4;
pTest_Cond = NaN(testNum, 2);
for i = 1:testNum
pTest_Cond(i, 1) = pScene_spatial_Cpost(XTest_Human_Cond(:, (i-1)*2 + 1), VTest_Human_Cond(:,(i-1)*2 + 1,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
pTest_Cond(i, 2) = pScene_spatial_Cpost(XTest_Human_Cond(:, (i-1)*2 + 2), VTest_Human_Cond(:,(i-1)*2 + 2,:), Xh20, Vh20, ZpostH20, RpostH20, YpostH20, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);
end