 function [Zpost, Ypost, Rpost, muCTpost, SigmaCTpost, runTimes, Kmax, wburn, wsample] = ...
    wood_ibp_learning_frontend(X,V,lda,eps,sigmaU,sigmaV,phi,sigmaC,alpha,bAlpha,bBeta, Kmax, wburn, wsample, stepNo, burnIn )
    rng(1);


    % Dependencies: wood_make_gibbs_y_spatial
    %               wood_make_gibbs_z_spatial

    initLatentNo = 1; %** 1 latent chunk at the beginning

    %** X is the binary observation matrix, N observed variables (N = 2 shapes), T trials (T = 10)
    [N, T] = size(X);

    %** Z[N x K]: binary dependency matrix, N observed variables (2 shapes), K latent variables (K chunks)
    %** Y[K x T]: binary latent variable state matrix, K latent variables (chunks), T trials (10)
    %** R[K]: latent variables' bias ~ beta(bAlpha, bBeta)

    Znew = zeros(N, initLatentNo);
    Rnew = betarnd(bAlpha, bBeta, 1, initLatentNo); 
    %Ynew = zeros(initLatentNo, T);
    Ynew = ((rand(initLatentNo, T) < repmat(Rnew', 1, T)) + 0);
    
    Zpost=cell(1,stepNo-burnIn); 
    Ypost=cell(1,stepNo-burnIn);
    Rpost=cell(1,stepNo-burnIn);
    muCTpost=cell(1,stepNo-burnIn);
    SigmaCTpost=cell(1,stepNo-burnIn);
    runTimes = NaN(stepNo, 1);

    for i = 1:stepNo
      tic
      display(i-1);  display(Znew); display(Rnew); %display(Ynew);

        %** with the updated states calculate the state
        % updating Y state matrix: update the states
        Ynew = wood_make_gibbs_y_spatial(Ynew,X,V,Znew,Rnew, ...
            lda, eps, sigmaU, sigmaV, phi, sigmaC);

        [Znew, Ynew, Rnew] = wood_make_gibbs_z_spatial(Ynew,X,V,Znew,Rnew, ...
            lda, eps, sigmaU, sigmaV, phi, sigmaC, bAlpha, bBeta, alpha, Kmax, wburn, wsample);

        [pi, muCTposti, SigmaCTpostInvi, SigmaCTposti] =  calc_pv_training(X, V, Ynew, Znew, sigmaU, sigmaV, phi, sigmaC, 0, 1);

        if (i>burnIn)
            
            Zpost{1,i-burnIn}=Znew;
            Ypost{1,i-burnIn}=Ynew;
            Rpost{1,i-burnIn}=Rnew;
            muCTpost{1,i-burnIn}=muCTposti;
            SigmaCTpost{1,i-burnIn}=SigmaCTposti;

        end
        time_i = toc
        runTimes(i) = time_i;
    end