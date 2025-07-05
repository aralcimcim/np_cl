function [p, muCTpost, SigmaCTpostInv, SigmaCTpost, zc, Psit] = ...
    calc_pv_training(x,v,y,Z, sigmaU, sigmaV, phi, sigmaC, varargin)
    %** calc_pv_training(X,V,Ynew,Z, sigmaU, sigmaV, phi, sigmaC,0,1);
    %** x: [N x T]      observation mtx
    %** v: [N x T x 2]  position of shapes
    %** y: [K x T]      state mtx of chunks
    %** Z: [N x K]      dependency mtx
    %** sigmaU: 1       variance of the chunk's position 
    %** sigmaV: 1       variance of shape's position, EQ 4 Orban_Wood
    %** sigmaC: 1       variance of relative position
    %** varagin:        1st: control of the dimensions to be evaluated (1: x only, 2: y only, 0 default: both)
    %**                 2nd: is logharithmic

    % Validated through: calc_pv_c
    % Validation example: p1=calc_pv_c(x(:,1),v(:,1,:),y(:,1),Z,c, 1, 2, 8); p2=calc_pv_c(x(:,2),v(:,2,:),y(:,2),Z,c, 1, 2, 8); p1*p2
    %
    % Dependencies: -> calc_pv_trial
    %
    % Optional arguments:
    %   dims - dimensions of the data to be evaluated: 
    %           1, x-dimension only
    %           2, y-dimension only
    %           0, both x&y dimensions (default)

    %** position variable of the observed feature, EQ 4 in Orban_Wood 2015

    args=varargin;
    nargs=length(args);
    if (nargs>0),
        dims=args{1};
    else
        dims=0;
    end
    if (nargs>1),
        logOutput=args{2};
    else
        logOutput=0;
    end

    [N T] = size(x);
    K = size(Z,2);
    %** Z[N x K] dependency, availLinkNo counts the dependencies between chunks and shapes
    %** note that here the dimensions are not N or K, but the number of all the
    %** dependecies
    availLinkNo = sum(Z(:));
    zc = false(N, K);
    zl = logical(Z);

    SigmaCtllhInv = zeros(availLinkNo, availLinkNo, T); %** EQ 44 Orban_Wood, [2 x 2 x 10]
    wt = zeros(availLinkNo, T, 2);                      %** [2 x 10 x 2]
    Psit = zeros(1, T, 2);                              %** Psi: EQ 54, def on p. 16, Orban_Wood
    for t = 1:T, %** for all trials calculate the inverse covariance mtx (EQ 44), Psit (EQ 54), and zct [N x K] logical
       [SigmaCtllhInv(:,:,t), wt(:,t,:), Psit(:,t,:), zct] = ...
            calc_pv_trial(x(:,t), v(:,t,:), y(:,t), Z, sigmaU, sigmaV, phi);
        zc = (zc|zct);                                  %[N, K], 0 if never active (x*y*z)
    end
    
    SigmaCTllhInv = sum(SigmaCtllhInv, 3);                                      %** sum for all trial EQ 48 O_W
    SigmaCTpostInv = SigmaCTllhInv + diag(ones(1, availLinkNo)) / (sigmaC ^ 2); %** EQ 50 in Orban_Wood

    SigmaCTpost = inv(SigmaCTpostInv);
    muCTpost = SigmaCTpost * sum(wt(:,:,1), 2);                                  %** EQ 51 O_W
    muCTpost(:,:,2) = SigmaCTpost * sum(wt(:,:,2), 2);

    muCTpostX = muCTpost(:,:,1);
    muCTpostY = muCTpost(:,:,2);
    
    if (sum(zc(:))), %** EQ 57
        nfx = mvnormpdfln( muCTpostX(zc(zl)), zeros(sum(zc(:)), 1), 'inv', SigmaCTpostInv(zc(zl), zc(zl)) );
        nfy = mvnormpdfln( muCTpostY(zc(zl)), zeros(sum(zc(:)), 1), 'inv', SigmaCTpostInv(zc(zl), zc(zl)) );
    else
        nfx = 0;
        nfy = 0;
    end

    %% ** calculating p(v_t | sigmaV) based on EQ 57 in O_W
    if (logOutput==0),
        p1=prod(Psit(:,:,1)) * (2*pi)^(-sum(zc(:))/2) * sigmaC^(-sum(zc(:))) ... 
            / exp(nfx);
        p2=prod(Psit(:,:,2)) * (2*pi)^(-sum(zc(:))/2) * sigmaC^(-sum(zc(:))) ...
            / exp(nfy);

        if (dims==0),
            p=p1*p2;
        elseif (dims==1),
            p=p1;
        elseif (dims==2),
            p=p2;
        else
            error('calc_pv_training: wrong setting for optional argument dims');
        end
    else
        lp1 = sum(log(Psit(:,:,1))) + ...
            log( (2*pi)^(-sum(zc(:))/2) * sigmaC^(-sum(zc(:))) ) - nfx ;
        lp2 = sum(log(Psit(:,:,2))) + ... 
            log( (2*pi)^(-sum(zc(:))/2) * sigmaC^(-sum(zc(:))) ) - nfy ;

        if (dims==0),
            p = lp1 + lp2;
        elseif (dims == 1),
            p = lp1;
        elseif (dims == 2),
            p = lp2;
        else
            error('calc_pv_training: wrong setting for optional argument dims');
        end

    end

end


