function p = pScene_ZR_spatial_Cpost(X, V, Z, R, Ypost, Xtraining, Vtraining, lda, eps, sigmaU, sigmaV, phi, sigmaC, Rtreshold, Vswitch)
    % X is a column vector e.g N=6, X=[1; 1; 1; 0 ; 1; 0; 0; 0; 1; 1; 0; 0]
    % this funtion calculates P(Scene | Z, R) = 
    % = sum_{all possible Y at trial T} P(X | Z, Y_T, R) * P(Y_T | R)
    
    if ~exist('Vswitch','var')
     % third parameter does not exist, so default it to something
      Vswitch = 1;
    end
    
    mapp = R > Rtreshold;
    Z = Z(:,mapp);
    R = R(:,mapp);
    Ypost = Ypost(mapp,:);
    
    K = size(Z, 2);
    [N, T] = size(X);
    [pv, CpostMu, CpostSigmaInv, CpostSigma, zc, Psit] = calc_pv_training(Xtraining, Vtraining, Ypost, Z, sigmaU, sigmaV, phi, sigmaC, 0,1);
    
    
    
%     infoMtx_1 = [Z V(:,:,1) V(:,:,2)];
%     infoMtx_2 = [CpostMu(:,:,1) CpostMu(:,:,2)];
%     display(infoMtx_1)
%     display(infoMtx_2)
  
    
    if K < 21 
        
        j = 0:2^K - 1;
        allY = dec2bin(j) - '0';
        lpXY = NaN(1, 2^K);
       
       %for i = 0:2^K-1
            Y = allY;%(i+1, :);
            lpV = NaN(size(Y,1),1);
            lpY =  sum(log(R.^Y') + log((1-R).^(1-Y)'));
            lpX = sum(log((1 - ((1-lda).^(Z*Y'))*(1-eps)) .* X + (((1-lda).^(Z*Y'))*(1-eps)) .* (1-X)));
                        
            for yi = 1:size(Y,1)
                %psi comes from the scnee, CPosts come from the training
               % [PsitScene, SigmaCTllhInv, muCTllh] = calc_pv_training_Scene(X, V, (Y(yi,:))', Z, sigmaU, sigmaV, phi, sigmaC);
                y = (Y(yi,:))';
                tmpPv =  calc_pv_trial_Cpost(X, V, y, Z, sigmaU, sigmaV, phi, CpostMu, CpostSigmaInv);
                if tmpPv == -Inf && (Vswitch==0)
                    lpV(yi) = 0;
                else lpV(yi) = tmpPv * Vswitch;
                end
                
                %calcPos_Cpost(X, V, y, Z, CpostMu, CpostSigmaInv, sigmaU, gridMin, gridMax);
            end
            lpXY = lpX + lpV' + lpY; 
       % end

      %  logP = 
        p = sum(exp(lpXY));  % here we can loose some precision
    
    else
        sizeMax = 2^K - 1;
        % first 20
        unit2 = 2^20;
        p_i = NaN(2^(K-20),1);
        
        %parfor n = 1:2^(K-20)
        for n = 1:2^(K-20)
            j = [((n-1) * unit2):(n*unit2 - 1) sizeMax];
            allY = dec2bin(j) - '0';
            %sizeJ = size(j,2);
            lpXY = NaN(1, size(j,2));
            
          %  for i = 0:size(j,2) - 2
                Y = allY;%(i+1, :);
                lpV = NaN(size(Y,1),1);
                lpY =  sum(log(R.^Y') + log((1-R).^(1-Y)'));
                lpX = sum(log((1 - ((1-lda).^(Z*Y'))*(1-eps)) .* X + (((1-lda).^(Z*Y'))*(1-eps)) .* (1-X)));
                for yi = 1:size(Y,1)
                    
                     lpV(yi) = calc_pv_trial_Cpost(X, V, (Y(yi,:))', Z, sigmaU, sigmaV, phi, CpostMu, CpostSigmaInv) * Vswitch;
                   % lpV(yi) = calc_pv_training(X, V, (Y(yi,:))', Z, sigmaU, sigmaV, phi, sigmaC,0,1);
                end
                lpXY = lpX + lpV' + lpY; 
          %  end
            p_i(n) = sum(exp(lpXY(1:end-1))); 
        end
        %  second: 21-25, hopefully it's enough, but i need to rewrite this
        
        
%         for n = (21+1):K
%             
%            
%             j = [2^(n-1):2^(n) - 1 sizeMax];
%             allY = dec2bin(j) - '0';
%             lpXY = NaN(1, size(j,2));
% 
%             for i = 0:size(j,2)-1
%                 Y = allY(i+1, :);
%                 lpY =  sum(sum(log(R.^Y) + log((1-R).^(1-Y))));
%                 lpX = sum(log((1 - ((1-lda).^(Z*Y'))*(1-eps)) .* X + (((1-lda).^(Z*Y'))*(1-eps)) .* (1-X)));
%                 lpXY(i+1) = lpX + lpY; 
%             end
%             pi(n-20) = sum(exp(lpXY)); 
% 
%            
%         
%         end
        
% for sizeMax:
    Y_sizeMax = dec2bin(sizeMax) - '0';
    lpY =  sum(sum(log(R.^Y_sizeMax) + log((1-R).^(1-Y_sizeMax))));
    lpX = sum(log((1 - ((1-lda).^(Z*Y_sizeMax'))*(1-eps)) .* X + (((1-lda).^(Z*Y_sizeMax'))*(1-eps)) .* (1-X)));
    
    lpV = calc_pv_trial_Cpost(X, V, Y_sizeMax', Z, sigmaU, sigmaV, phi, CpostMu, CpostSigmaInv) * Vswitch;
    %lpV = calc_pv_training(X, V, Y_sizeMax', Z, sigmaU, sigmaV, phi, sigmaC,0,1);

    lpXY= lpX + lpV + lpY; 

    p_sm = sum(exp(lpXY));
    p = sum(p_i) + p_sm;
        
        
        
    end
end