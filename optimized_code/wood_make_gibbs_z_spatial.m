function [Znew, Ynew, Rnew] = wood_make_gibbs_z_spatial(Y,X,V,Z,R, lda, eps, sigmaU, sigmaV, phi, sigmaC, Balpha, Bbeta, alpha,Kmax, wburn, wsample)
%rng(0);
% Dependencies: -> calc_pv_training

[N, K] = size(Z);
[T] = size(X, 2);
 
Znew = Z;
lpv_y = zeros(1, 2); 
lpx_y = zeros(1, 2);

wall = wburn + wsample; 

for n = 1:N %** for each shapes
    
    %% update Z_n for all K
    Znew = Z;
    lpv_y_new = calc_pv_training(X, V, Y, Znew, sigmaU, sigmaV, phi, sigmaC, 0, 1);  %** p(v | sigma_v) EQ 57 OW
    lpx_y_new = calc_px_training(X(n,:), Y, Znew(n,:), lda, eps, 1);                 %** p(x | y, Z) EQ 3 OW
     
    for k = 1:K %** for each chunks
        if (n < N) %** m_{-n,k}
            m_mnk = sum(Znew([1:n-1 n+1:N], k)); 
        else %** if n=N, i.e. last shape 
            m_mnk = sum(Znew(1:n-1, k)); 
        end
        
        %**! comment m_mnk
        
        th_k = m_mnk/N; %** P(z_ik | Z_{-i, k}) p.5 OW

        if (m_mnk > 0)
            %** z_nk = 1/0 based on Z (in EQ 12 OW)
            lpv_y(1, 1+Znew(n,k)) = lpv_y_new;%calc_pv_training(X, V, Y, ZAct, sigmaU, sigmaV, phi, sigmaC,0,1); %lpv_y_new;
            lpx_y(1, 1+Znew(n,k)) = lpx_y_new;%calc_px_training(X(n,:), Y, ZAct(n,:), lda, eps, 1); 
           
            %** z_nk = 1/0, the opposite of Z (in EQ 12 OW)
            Znew(n, k) = 1 - Znew(n, k);
            lpv_y(1, 1 + Znew(n,k)) = calc_pv_training(X,V,Y,Znew, sigmaU, sigmaV, phi, sigmaC,0,1); 
            lpx_y(1, 1 + Znew(n,k)) = calc_px_training(X(n,:),Y, Znew(n,:), lda, eps, 1);
            
            %** EQ 12 OW
            lPznk1_XZmnkY = log(th_k) + lpx_y(1, 2) + lpv_y(1, 2); 
            lPznk0_XZmnkY = log(1 - th_k) + lpx_y(1, 1) + lpv_y(1, 1);
            %** EQ 13 OW
            rt = exp(lPznk1_XZmnkY - lPznk0_XZmnkY);
            
            %Pznk_XZmnkY = p(Znk = 1)/(p(Znk = 1) + p(Znk = 0))
            if (rt == Inf)
                Pznk_XZmnkY = 1;
            else
                Pznk_XZmnkY = rt/(rt+1); %** similarly to y_spatial, this is P1/(P1+P0) in EQ 13 (without normalizing the probabilities), EQ 13 OW
            end
            
        else
            Pznk_XZmnkY = 0;
        end
        
        uni = rand(1);
        Znew(n, k) = Pznk_XZmnkY > uni; %** Znew(shape,chunk) is a random binary with Bern ~ Pznk_XZmnkY=p in EQ 13 OW
        
        lpv_y_new = lpv_y(1, Znew(n,k) + 1);
        lpx_y_new = lpx_y(1, Znew(n,k) + 1);
    
    end
    
    Z = Znew;
    
    %% update R for all K
    for k = 1:K
        R(1,k) = betarnd(Balpha + sum(Y(k,:)), Bbeta + T - sum(Y(k,:)));  %** posterior p(r | T, #active) ~ Beta(a + #active, b + T - #active) 
    end
       
    %% Heading for possible new latents
   
    PKnew_XZY = zeros(1,Kmax+1);
    PKnew_XZY_an = zeros(1,Kmax+1);
  
    for Knew = 0:Kmax %** running the scenario for 0, 1, .. Kmax new chunks
        lpy_xz = zeros(1,wsample);
        ZAct = [Z zeros(N,Knew)]; ZAct(n, K+1:K+Knew)=1;
        RAdd = betarnd(Balpha, Bbeta, 1, Knew);
        RAct = [R RAdd];
        YAct = [Y; ((rand(Knew, T) < repmat(RAdd', 1, T)) + 0)]; 
        
        lpv_y_act = zeros(1,2); 
        lpx_y_act = zeros(1,2);
        
        lpv_y_new = calc_pv_training(X, V, YAct, ZAct, sigmaU, sigmaV, phi, sigmaC, 0, 1);
        lpx_y_new = calc_px_training(X(n,:), YAct, ZAct(n,:), lda, eps, 1);
             
   %% this comes back now for vectorized p   
   
        if (Knew > 0)
            for wi = 1:wburn+wsample %** sampling 11 times, counting only the last one
                for t = 1:T %** for all trials
                    for k = 1:Knew %** calculate Y state values for the new chunks, same as in wood_make_gibbs_y_spatial
                        
                        lpv_y_act(1, 1 + YAct(K+k, t)) =  lpv_y_new;
                        lpx_y_act(1, 1 + YAct(K+k, t)) =  lpx_y_new;
                        
                        YAct(K+k, t) = 1 - YAct(K+k, t);
                        lpv_y_act(1, 1 + YAct(K+k, t)) = calc_pv_training(X,V,YAct,ZAct, sigmaU, sigmaV, phi, sigmaC,0,1);%...
                        lpx_y_act(1, 1 + YAct(K+k, t)) = calc_px_training(X(n,:), YAct, ZAct(n,:), lda, eps, 1);  %** CHANGED from calc_px_training(X,YAct,ZAct, lda, eps, 1);
                        
                        p = RAct(1, K + k);
                        %YAct(K+k,t)=1;
                        lPYk1_ZXYmk = log(p) + lpx_y_act(1,2) + lpv_y_act(1,2);
                        %YAct(K+k,t)=0;
                        lPYk0_ZXYmk = log(1-p) + lpx_y_act(1,1) + lpv_y_act(1,1);

                        rt = exp(lPYk1_ZXYmk - lPYk0_ZXYmk );
                        
                        if (rt == Inf) 
                            PYk1_ZXYmk = 1;
                        else 
                            PYk1_ZXYmk = rt / (rt+1); 
                        end
                        
                        uni = rand(1,1);
                        YAct(K+k,t) = ...
                            ( PYk1_ZXYmk > uni );

                        lpv_y_new = lpv_y_act(1, YAct(K + k, t) + 1);
                        lpx_y_new = lpx_y_act(1, YAct(K + k, t) + 1);
                        
                    end
                end
                 
                for k = 1:Knew %** calculate R bias values for the new chunks, EQ
                    RAct(1, K + k) = betarnd(Balpha + sum(YAct(K + k, :)), Bbeta + T - sum(YAct(K + k, :))); %** posterior p(r | T, #active) ~ Beta(a + #active, b + T - #active) 
                end
               
                              
                if (wi>wburn) %** in the sample steps
                                                          
                    lpy_v(wi-wburn) = lpv_y_new;
                    lpy_x(wi-wburn) = lpx_y_new;   
                    lpy_xv(wi-wburn) = lpx_y_new + lpv_y_new ;
                end
                  
            end        
            
        else
            lpx_y_new = calc_px_training(X(n,:), YAct, ZAct(n,:), lda, eps, 1);
            lpv_y_new = calc_pv_training(X,V,YAct,ZAct, sigmaU, sigmaV, phi, sigmaC,0,1);
            
            lpy_xv = lpv_y_new + lpx_y_new ;                   
                      
        end
        
        % LOG harmonic mean
        % ** !! comment on importance sampling
        if Knew == 0 
            py_xzKnew_5 = lpy_xv; 
        else
            py_xzKnew_5 = log(wsample) - log(1) + max(lpy_xv) - log(sum(1./exp(lpy_xv - max(lpy_xv))));
        end
       
        PKnew_XZY_5(Knew + 1) = (py_xzKnew_5 + log(poisspdf(Knew, alpha/N)));
        Ysample{Knew + 1} = YAct;
        Rsample{Knew + 1} = RAct;
    end
        
    lpdf_5 = PKnew_XZY_5 - max(PKnew_XZY_5);
    pdf_5 = exp(lpdf_5)/sum(exp(lpdf_5));

    bar = rand(1);
    j = 1; while (sum(pdf_5(1:j))<bar), j = j+1; end
    Knew = j - 1; %** here we've decided on the new chunk number for shape n
    
    Znew = [Z zeros(N, Knew)]; 
    Znew(n, K+1:K+Knew) = 1;
    
    %display(n);
    Z = Znew;
    Y = Ysample{Knew + 1}; 
    R = Rsample{Knew + 1};
        
    nodes2keep = find(sum(Z,1));
    if(~isempty(nodes2keep) )
        Z = Z(:, nodes2keep);
        Y = Y(nodes2keep, :);
        R = R(1,nodes2keep, :);
        K = size(Z, 2);
    elseif isempty(Z)
        Z = zeros(N,1);
        Y = zeros(1,T);
        R = 0;
        K = size(Z,2);  
    else    
        Z = zeros(N,1);
        Y = zeros(1,T);
        R = 0;
        K = size(Z,2);
    end       
         
end

Znew = Z;
Ynew = Y;
Rnew = R;