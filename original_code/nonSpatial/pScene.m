function p = pScene(XScene, Zpost, Rpost, lda, eps)
    % X is a column vector
    % Zpost and Rpost are the cell arrays, outputs of the IBP algorithm
    % lda and eps are the noisy-or (input) parameters of the IBP model
    
    % pScene calculates the probability P(XScene) =
    % = sum_{all output datum: for all steps} P(XScene | Z, R) * P(Z, R) =
    % = sum pScene_ZR(XScene, Z, R..) * 1/stepNo
    
    stepNo = size(Zpost, 2);
    pCond = NaN(stepNo,1);
    
    for i = 1:stepNo
        Z = Zpost{i};
        R = Rpost{i};
        pCond(i) = pScene_ZR(XScene, Z, R, lda, eps);        
    end
    
    p = sum(pCond) * 1/stepNo;

end