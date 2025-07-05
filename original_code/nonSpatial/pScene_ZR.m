function p = pScene_ZR(X, Z, R, lda, eps)
    % X is a column vector e.g N=6, X=[1; 1; 1; 0 ; 1; 0; 0; 0; 1; 1; 0; 0]
    % this funtion calculates P(Scene | Z, R) = 
    % = sum_{all possible Y at trial T} P(X | Z, Y_T, R) * P(Y_T | R)
    
    K = size(Z, 2);
    
    if K < 21 
        
        j = 0:2^K - 1;
        allY = dec2bin(j) - '0';
        lpXY = NaN(1, 2^K);

       
       %for i = 0:2^K-1
            Y = allY;%(i+1, :);
            lpY =  sum(log(R.^Y') + log((1-R).^(1-Y)'));
            lpX = sum(log((1 - ((1-lda).^(Z*Y'))*(1-eps)) .* X + (((1-lda).^(Z*Y'))*(1-eps)) .* (1-X)));
            lpXY = lpX + lpY; 
       % end

      %  logP = 
        p = sum(exp(lpXY));  % here we can loose some precision
    
    else
        sizeMax = 2^K - 1;
        % first 20
        unit2 = 2^20;
        p_i = NaN(2^(K-20),1);
        
        parfor n = 1:2^(K-20)
        
            j = [((n-1) * unit2):(n*unit2 - 1) sizeMax];
            allY = dec2bin(j) - '0';
            %sizeJ = size(j,2);
            lpXY = NaN(1, size(j,2));

          %  for i = 0:size(j,2) - 2
                Y = allY;%(i+1, :);
                lpY =  sum(log(R.^Y') + log((1-R).^(1-Y)'));
                lpX = sum(log((1 - ((1-lda).^(Z*Y'))*(1-eps)) .* X + (((1-lda).^(Z*Y'))*(1-eps)) .* (1-X)));
                lpXY = lpX + lpY; 
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
lpXY= lpX + lpY; 
       
p_sm = sum(exp(lpXY));
p = sum(p_i) + p_sm;
        
        
        
end