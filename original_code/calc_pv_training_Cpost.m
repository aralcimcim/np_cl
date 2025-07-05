function [p, muCTpost, SigmaCTpostInv, SigmaCTpost] = ...
    calc_pv_training_Cpost(x,v,y,Z, muCTpost, SigmaCTpostInv, zc, Psit,  SigmaCTllhInv, muCTllh, sigmaU, sigmaV, phi, sigmaC, varargin)
   
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
   
    zl = logical(Z);
   
    muCTllhX = muCTllh(:,:,1);
    muCTllhY = muCTllh(:,:,2);
    
    if (sum(SigmaCTllhInv(:) == 0)) 
         SigmaSumInv = inv(inv(SigmaCTpostInv) + SigmaCTllhInv);
    else SigmaSumInv = inv(inv(SigmaCTpostInv) + inv(SigmaCTllhInv));
    end
    
    
    if (sum(zc(:))), %** EQ 57
        nfx = mvnormpdfln( muCTllhX(zc(zl)), muCTpost(zc(zl)), 'inv', SigmaSumInv(zc(zl), zc(zl)) );
        nfy = mvnormpdfln( muCTllhY(zc(zl)), muCTpost(zc(zl)), 'inv', SigmaSumInv(zc(zl), zc(zl)) );
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


