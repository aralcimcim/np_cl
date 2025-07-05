function [Z,K] = ibp_generate(N,alpha,K_specified)
% function [Z,K] = ibp_generate(N,alpha,K_specified)
%
%   WARNING!  This function is not for the faint of heart!  It produces a
%   matrix Z according to the Indian Buffet Process (with the dirichlet
%   hyperparameter alpha) by generated a bunch of random matrixes until 
%   one of the correct dimensionality K_specified pops up.  if K_specified
%   isn't specified then the function returns a matrix of output
%   dimensionality N but random K.  Defaults: alpha = 2.2 K_specified = -1;
%   (which also means unspecified)
switch(nargin)
    case 0
        N = 20;
        K_specified = -1;
        alpha = 2.2;
    case 1
        K_specified = -1;
        alpha = 2.2;
    case 2
        K_specified = -1;
end

debug = 0;

tempK = 100;
K = Inf;


while(K ~= K_specified)

%     rand('state',0);

    Z = zeros(N,tempK);


    nfirst = poissrnd(alpha,1);
    if(nfirst == 0)
        continue;
    end
    Z(1,1:nfirst) = 1;
    m = sum(Z);
    max_c = nfirst;

    for(i=2:N)
        rns = rand(1,tempK);
        pr = m/i;
        v = pr>rns;
        Z(i,:) = v;
        %       Z(i,1:max_c) = v(1:max_c);

        k_more = poissrnd(alpha/i,1);
        if(k_more >=1)
            Z(i,max_c:(max_c+k_more)) = 1;
        end
        max_c = max_c+k_more;
        m = sum(Z);
    end

    K = max_c;
    if(debug)
        disp(['K = ' num2str(K)]);
    end
    Z = Z(:,1:K);

    if(K_specified==-1)
        break;
    end
end