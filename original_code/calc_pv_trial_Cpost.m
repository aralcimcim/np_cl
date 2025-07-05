function [p, SigmaCtllhInv, wt, Psit, zc] = calc_pv_trial(x, v, y, Z, sigmaU, sigmaV, phi, CpostMu, CpostSigmaInv)
    
   
   [N, K] = size(Z);
   actShape = sum(x);
   ZAct = Z.*repmat(y',N,1);               
   zc = logical(ZAct.*repmat(x,1,K));           
   zl = logical(Z);
   
   
    ZAct = Z.*repmat(y',N,1);                   %** ZAct [N x K] is the dependency mtx for the active chunks (in the trial) only
    phi = phi.*ZAct;                            %** phi: [N x K],  0 if no dependency or chunk state is 0 in the trial,  phi if the shape is dependent and chunk is present in the trial
    %(phi + rand(N,K)).*ZAct; 
    zc = logical(phi.*repmat(x,1,K));           %** zc:  [N x K],   0 if shape is not present, or no dependency or chunk state is 0 in the trial, 1 if shape is present, there is dependency, and chunk state is 1 in the trial
    
    xAct = find(x);                             %** indices of 1 in the x[N x 1] vector, i.e. which shapes are present in the given trial
    xActNo = length(xAct);                      %** how many shapes are present in the trial
    A = phi ./ (1 + repmat(sum(phi,2),1,K));    %** A [N x K]: 0, if no dependency, or chunk state is 0 in the trial, Aij (e.g EQ 21) otherwise
    sSq = sigmaV^2./(1 + sum(phi,2));           %** sSq = 1 / (1 + 4 * #linked active chunk)    if the shape is dependent and chunk is present in the trial
    
    
    ptAct = find( sum(ZAct.*repmat(x, 1, K), 1) );  %** ptAct is the serial num of the chunks with active links to shapes in the trial
    ptNo = length(ptAct);                           %** how many active chunks in the trial 
    SigmaInv = zeros(ptNo, ptNo, length(xAct));
        
    omegaBU = zeros(N);                             %** expressed in EQ 36 O_W
    omega = 1 / sigmaV^4 * (diag(sSq(xAct)));  
                                                       

      %** CHANGE: put this 4 lines into an if to avoid errors  
     availLinks=[];  availLinksRow=[]; availLinksCol=[];
     if sum(sum(Z)) ~= 0
        [availLinksRow(:,1) availLinksCol(:,1)] =  find(Z); %** position of the nonzero links, first row (1,1) means that the 1st shape has link to the 1st chunk
        availLinks = [availLinksRow availLinksCol];
        availLinkNo = size(availLinks, 1); 
        SigmaCtllhInv = zeros(availLinkNo);
        wt = zeros(availLinkNo,1,2); 
     else availLinkNo = 0;
        SigmaCtllhInv = 0; %** !!! this is probably wrong just to make the code run
        wt = 0;             %** !!! this is probably wrong just to make the code run
     end
     %** end of CHANGE


    if (ptNo > 0),    

        for xi = 1:length(xAct),  %** for all shapes present in the trial
            SigmaInv(:,:,xi) = ( A(xAct(xi), ptAct)' * A(xAct(xi), ptAct) ) / sSq(xAct(xi)); %** O_W p.10 sigma_u_j_i'_k'^-1
        end
        
        SigmaCVInv = sum(SigmaInv, 3);                         %** O_W p.8
        SigmaCInv = SigmaCVInv + 1 / (sigmaU^2) * eye(ptNo);   %** O_W p.10 sigma_c^-1
        SigmaC = inv(SigmaCInv);

        omega = omega - 1/sigmaV^4 * ((A(xAct,ptAct) * (SigmaC * A(xAct,ptAct)'))); %**[N x K]: EQ 36 in O_W
        omegaBU = zeros(N);
        omegaBU(xAct, xAct) = omega;
        
        for i = 1:availLinkNo, %** O_W EQ 44
            for j = 1:availLinkNo,
                if (x(availLinks(i,1)) == 1 && x(availLinks(j,1)) == 1), %** if there is link to the shape and the shape is presented in the trial
                    SigmaCtllhInv(i,j) = ...
                        omegaBU(availLinks(i,1), availLinks(j,1)) * ...
                        phi(availLinks(i,1), availLinks(i,2)) * ...       
                        phi(availLinks(j,1), availLinks(j,2));
                end
            end

            vnan = v;
            vnan(isnan(vnan)) = 0;
            wt(i, 1, 1) =   phi(availLinks(i, 1), availLinks(i, 2)) * ...  %** EQ 46 in O_W, handle nan values
                            omegaBU(availLinks(i, 1), :) * ...
                            ((1 + sum(phi, 2)) .* vnan(:, :, 1));
            wt(i, 1, 2) =   phi(availLinks(i, 1), availLinks(i, 2)) * ...
                            omegaBU(availLinks(i, 1), :) * ...
                            ((1 + sum(phi, 2)) .* vnan(:, :, 2));
        end
    else
       SigmaC = 1;
    end
      
    Psit = zeros(1,1,2);                    %** EQ 54 O_W, to simplify the calculation of P(v_t | C, sigmaV)

    vp = ((1 + sum(phi, 2)) .* v(:,:,1));   %** 1st coordinate of shape in trial * (1 + 4 * #linked active chunk)
    vpnan = vp;                             %** added by me, handling nans in vp
    vpnan(isnan(vpnan)) = 0;
    Psit(:,:,1) =   (2*pi)^(-xActNo/2) * sigmaU^(-ptNo) * ...           %** in EQ 54 
                    1/sqrt(prod(sSq(xAct))) * det(SigmaC)^(.5) * ...    %** sqrt because sj^2 is calculated in the form of sSq
                        exp(-.5 * (vpnan' * (omegaBU * vpnan)));

    vp = ((1 + sum(phi,2)) .* v(:,:,2));
    vpnan = vp;
    vpnan(isnan(vpnan)) = 0;
    Psit(:,:,2) =  (2*pi)^(-xActNo/2) * sigmaU^(-ptNo) * ...
                    1/sqrt(prod(sSq(xAct))) * det(SigmaC)^(.5) * ...
                    exp(-.5 * (vpnan' * (omegaBU * vpnan)));
    
    %% addendum for the actual dimensions in the trial
    dim = sum(zc(zl));
    SigmaCtllhInv_trial = SigmaCtllhInv(zc(zl), zc(zl));
    wt_trial = wt(zc(zl), :, :);
    
    
    CpostMu_trial = CpostMu(zc(zl), :, :);
    CpostSigmaInv_trial = CpostSigmaInv(zc(zl), zc(zl));
    
    if det(SigmaCtllhInv_trial)<1e-15
        SigmaCtllhInv_trial = SigmaCtllhInv_trial + eye(dim)/100000;
    end
    muCtllh_trial(:,:,1) = inv(SigmaCtllhInv_trial) * wt_trial(:,:,1);
    if isempty(wt_trial) muCtllh_trial(:,:,2) = muCtllh_trial(:,:,1); %% this is a technical adjustment for handling empty stuff
    else muCtllh_trial(:,:,2) = inv(SigmaCtllhInv_trial) * wt_trial(:,:,2); 
    end
    
    sigmaAggrInv = inv(inv(SigmaCtllhInv_trial) + inv(CpostSigmaInv_trial));
          
    lnorm_x = mvnormpdfln( muCtllh_trial(:,:,1), CpostMu_trial(:,:,1), 'inv', sigmaAggrInv);
    lnorm_y = mvnormpdfln( muCtllh_trial(:,:,2), CpostMu_trial(:,:,2), 'inv', sigmaAggrInv);
    
    pv_x = log(Psit(:,:,1)) + log((2*pi)^(dim/2)) + log(det(inv(CpostSigmaInv_trial))^0.5) + (0.5 * muCtllh_trial(:,:,1)' * SigmaCtllhInv_trial * muCtllh_trial(:,:,1)) + lnorm_x;
    pv_y = log(Psit(:,:,2)) + log((2*pi)^(dim/2)) + log(det(inv(CpostSigmaInv_trial))^0.5) + (0.5 * muCtllh_trial(:,:,2)' * SigmaCtllhInv_trial * muCtllh_trial(:,:,2)) + lnorm_y;
                
                
    p = pv_x + pv_y;            
                
    
end

