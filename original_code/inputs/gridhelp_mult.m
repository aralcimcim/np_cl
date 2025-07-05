% function for plotting the first T scenes from the input
%
% X is the shape matrix, in case of 12 shapes and 144 trials the
% dimensions are 12x144
% V is the position matrix, dimensions are 12x144x2
%
% dependencies: lscatter.m 
% example usage: gridhelp_mult(trainingdata.shapes, trainingdata.pos, 144/4)

function gridhelp_mult(X,V,scenes)

    plotsInRow = 8;
    minV = -2;
    maxV = 2;

    numfigs = length(scenes);
    numcol = min(plotsInRow, numfigs);
    numrow = ceil(numfigs/plotsInRow);

    figure
    for t=1:numfigs
        sc=scenes(t);
        subplot(numrow, numcol, t);
        
        minV=min(V(:));
        maxV = max(V(:));
        N = size(X,1);

        lscatter(V(:,sc,1), V(:,sc,2), [1:N], 'FontSize', 10)
        axis([minV-0.5 maxV+0.5 minV-0.5 maxV+0.5])
        hold on

        grid0 = minV-0.5:1:maxV+0.5;
        grid1 = [grid0; grid0];
        grid2 = repmat([minV-0.5; maxV+0.5],1,length(grid0));
        plot(grid1,grid2,'k')
        plot(grid2,grid1,'k')

        tickValues2 = minV:1:maxV;
        
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
%         set(gca,'XTick',tickValues2);
%         set(gca,'YTick',tickValues2);
        title(strcat(num2str(sc)));
    end
    
end