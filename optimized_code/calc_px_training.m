function p = calc_px_training(X, Y, Z, lda, eps, varargin)
    %** X: [N x 1] = [12 x 1]           binary observation mtx for all shapes through the given trial, X = X(:,t)
    %** Y: [K x 1] = [K x 1]            binary state mtx for all chunks for the given trial, Y = Y(:,t)
    %** Z: [N x K] = [144 x K]          binary dependency mtx for all shapes and chunks 

    %** pp = p(x_n^t | y^t, Z), shape_n's state probability given the chunks and the Z dependencies
    %** state variable of the observed feature, EQ 3 in Orban_Wood 2015
    %** noisy-or: https://view.officeapps.live.com/op/view.aspx?src=http://pages.cs.wisc.edu/~dpage/cs731/lecture4.ppt

    %%
    args=varargin;
    nargs=length(args);
    if (nargs>0)
        isLog=args{1};
    else
        isLog=0;
    end

    %** pp is calculated for all shapes for the given trial
    pp = (1 - ((1-lda).^(Z*Y))*(1-eps)) .* X + ... %** Z * Y: [N x K] * [K x 1] binary, number of active chunks affecting the shape
        (((1-lda).^(Z*Y))*(1-eps)) .* (1-X);       %** the same calculation for the non-present shapes 

    %** sum(pp) represents the joint probability of the shapes state through the trial
    if (isLog==0)
        p=prod(pp(:));
    else
        pp=log(pp);
        p=sum(pp(:));
    end

end