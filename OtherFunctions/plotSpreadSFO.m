function outVals = plotSpreadSFO(inVals)

numCols = length(inVals);

allVals = [];
for colNum = 1:numCols
   allVals = cat(1,allVals,[colNum.*ones(length(inVals{colNum}),1),inVals{colNum}]);
   for barNum = 1:max(inVals{colNum}(:,1))
       barVals(colNum) = mean(inVals{colNum});
       errorBarVals(colNum) = std(inVals{colNum})./sqrt(length(inVals{colNum}));
   end
end

plotFig = figure('color','w');
plotAx = axes('parent',plotFig,'parent',plotFig);
hold(plotAx,'on');

bar(1:numCols,barVals,'facecolor','none','edgecolor','k','parent',plotAx);
errorbar(1:numCols,barVals,errorBarVals,'linestyle','none','color','k','marker','none',...
   'parent',plotAx);
plotSpread(allVals(:,2),'distributionIdx',allVals(:,1),'distributionMarkers','o')

% set(plotAx,'tickdir','out','xticklabel',{'3 mW','15 mW','15 mW + Cs'},...
%     'ylim',[-5 40],'ytick',0:10:40,'xlim',[0.25 numCols + 0.75]);
% ylabel(plotAx,'Current (pA)');
set(plotAx,'tickdir','out','xlim',[0.25 numCols + 0.75]);

[p,tbl,stats,terms] = anovan(allVals(:,2)',{allVals(:,1)},'model','interaction','varnames',{'recType'});

disp(p);