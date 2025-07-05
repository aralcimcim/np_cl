function Ynew = wood_make_gibbs_y_spatial(Y,X,V,Z,R, lda, eps, sigmaU, sigmaV, phi, sigmaC)
    %rng(0);

    % Dependencies: -> calc_pv_training

    [K, T] = size(Y); %** K chunks, T = 144 trials
    Ynew = Y; 

    %** position variable of the latent chunk
    lpv_y = zeros(1,2);
    %** state variable of the latent chunk
    lpx_y = zeros(1,2);

    %** note that p(V,params) depends on C over time, therefore we can't
    %break it down to time steps, we have to consider all t=1..T steps
    lpv_y_new = calc_pv_training(X, V, Ynew, Z, sigmaU, sigmaV, phi, sigmaC,0,1);
    lpx_y_new = calc_px_training(X(:,1), Ynew(:,1), Z, lda, eps, 1);

    %** we update the state for each chunks, for all trials.
    for t = 1:T, %** for all trials
        
        for k = 1:K, %** for all chunks

            %** changing the coordinates: 2nd if the chunk is relevant for the trial, 1st if not
            lpv_y(1, 1 + Ynew(k,t)) = lpv_y_new;
            lpx_y(1, 1 + Ynew(k,t)) = lpx_y_new;

            Ynew(k,t) = 1 - Ynew(k,t); %** flipping Ynew 0/1 value to change the remaining coordinates %** EQ 12 Woods 2006, this represents (1-a), but I don't get the usage 

            lpv_y(1, 1 + Ynew(k,t)) = calc_pv_training(X, V, Ynew, Z, sigmaU, sigmaV, phi, sigmaC,0,1);  %calc_pv_training(X,V,Ynew,Z, sigmaU, sigmaV, phi, sigmaC,0,1); 
            lpx_y(1, 1 + Ynew(k,t)) = calc_px_training(X(:,t), Ynew(:,t), Z, lda, eps, 1);

            % bias of the kth chunk ~beta(bAlpha, bBeta)
             p = R(1,k); %

            %** p(y = 1 | ...)
            lPYk1_ZXYmk = log(p) + lpx_y(1, 2) + lpv_y(1, 2);

            %** p(y = 0 | ...)
            lPYk0_ZXYmk = log(1 - p) + lpx_y(1, 1) + lpv_y(1, 1);

            rt = exp(lPYk1_ZXYmk - lPYk0_ZXYmk); %EQ 10 Orban_Woods %** EQ 12 Woods 2006, don't get the substraction, why not exp(lPYk1_ZXYmk + lPYk0_ZXYmk)

            %** Ynew(k, t) ~ Bern(rt/(rt+1)), this is because p1 + p0 should be
            %** normalized to sum up to 1, but instead we calculate p(y=1) as Bern~ p1/(p1+p0)
            %** QQ there should be a check to sum this up to 1, and if it fails than we have a problem
            if (rt == Inf) 
                PYk_ZXYmk = 1;
            else PYk_ZXYmk = rt / (rt + 1); 
            end;
           
            Ynew(k,t) = (PYk_ZXYmk > rand(1,1)); 
            
            lpv_y_new = lpv_y(1, Ynew(k, t) + 1);
            lpx_y_new = lpx_y(1, Ynew(k, t) + 1);
           
        end
    end
    debug=0;
end
