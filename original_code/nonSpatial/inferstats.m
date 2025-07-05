function [Ek, EZZt, in_degree_error, structure_error, in_degree_error_ratio, structure_error_ratio] = inferstats(Zsamples,Z,start,graphics,knum)
%** added knum input for plotting
if(nargin<3)
    start = 0;
end


if(nargin<4)
    graphics=0;
end

% if(~iscell(Zsamples))
%     foo = Zsamples;
%     Zsamples = cell(1);
%     Zsamples{1} = foo;
% end

ksum = 0;
ZZt = zeros(size(Zsamples{1}*Zsamples{1}'));
num_samples = length(Zsamples);

for(i=(start+1):num_samples)
    ksum = ksum+size(Zsamples{i},2);
    ZZt = ZZt + Zsamples{i}*Zsamples{i}';
end
Ek = ksum/(num_samples-start);
EZZt = ZZt/(num_samples-start);

if(nargin>1)
    truth = Z*Z';

    in_degree_error = sum(abs(diag(truth)-diag(EZZt)));
    structure_error = sum(sum(abs(triu(truth,1)-triu(EZZt,1))));
    
    in_degree_error_ratio = in_degree_error/sum(diag(truth));
    total_link_measure =(sum(sum(triu(truth,1))));
    if(total_link_measure~=0)
        structure_error_ratio = structure_error/total_link_measure;
    else
        structure_error_ratio = structure_error;
    end
    
    if(graphics)
        %** CHANGE from figure(17);
        figure(knum);
        colormap hot
        subplot(1,2,1)
        imagesc(Z*Z')
        title(strcat('ZZ'' Actual, K = ',num2str(knum))) %** added knum
        colorbar
        subplot(1,2,2)
        imagesc(EZZt)
        colorbar
        title(strcat('ZZ'' Estimated, K = ', num2str(knum))) %** added knum
    end
else
    in_degree_error=-1;
    structure_error=-1;
    in_degree_error_ratio=0;
    structure_error_ratio=0;
end