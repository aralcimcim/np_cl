function p = pScene_spatial_Cpost(XScene, VScene, Xtraining, Vtraining, Zpost, Rpost, Ypost, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch)
    % X is a column vector
    % Zpost and Rpost are the cell arrays, outputs of the IBP algorithm
    % lda and eps are the noisy-or (input) parameters of the IBP model
    
    % pScene calculates the probability P(XScene) =
    % = sum_{all output datum: for all steps} P(XScene | Z, R) * P(Z, R) =
    % = sum pScene_ZR(XScene, Z, R..) * 1/stepNo
    if ~exist('Vswitch','var')
     % third parameter does not exist, so default it to something
      Vswitch = 1;
    end
    
    stepNo = size(Zpost, 2);
    pCond = NaN(stepNo,1);
    
    for i = 1:stepNo
        Z = Zpost{i};
        R = Rpost{i};
        Y = Ypost{i};
%         display(i)
        pCond(i) = pScene_ZR_spatial_Cpost(XScene, VScene, Z, R, Y, Xtraining, Vtraining, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch);        
    end
  
    p = sum(pCond) * 1/stepNo;

end