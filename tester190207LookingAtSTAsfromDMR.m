fnames = what;
fnames = fnames.mat;

numFiles = length(fnames);
pvCounter = 1;
msnCounter = 1;
pvStore = [];
msnStore = [];
pvSpikeNum = [];
msnSpikeNum = [];
pvTuneStore = [];
msnTuneStore = [];
pvHistStore = [];
msnHistStore = [];

for i= 1:numFiles
    disp(i)
    load(fnames{i})
    numUnits = length(s.DesignationName);
    %find cell type target
    [indCellType] = functionCellStringFind(masterHeader,'CellType');
    [indSpikes] = functionCellStringFind(masterHeader,'NumSpikes');
    for j = 1:numUnits;
        %determine cell type
        if masterData(j,indCellType) == 1 %THIS IS FSI
            pvStore(:,:,pvCounter) = s.(s.DesignationName{j}).DMRAverage;
            pvSpikeNum(pvCounter) = masterData(j,indSpikes);
            pvTuneStore(:,:,pvCounter) = s.(s.DesignationName{j}).BinTone;
            pvHistStore(:,pvCounter) = s.(s.DesignationName{j}).AllHistograms;
            pvCounter = pvCounter + 1;
        elseif masterData(j,indCellType) == 0 %THIS IS MSN
            msnStore(:,:,msnCounter) = s.(s.DesignationName{j}).DMRAverage;
            msnSpikeNum(msnCounter) = masterData(j,indSpikes);
            msnTuneStore(:,:,msnCounter) = s.(s.DesignationName{j}).BinTone;
            msnHistStore(:,msnCounter) = s.(s.DesignationName{j}).AllHistograms;
            msnCounter = msnCounter + 1;
            
        end
    end
end


ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
load(ttlInfo)
%restrict data to roughly 4-32 kHz
pvStore = pvStore(19:36,:,:);
msnStore = msnStore(19:36,:,:);

%now lets try plotting things?


numMSNs = size(msnStore,3);
numPVs = size(pvStore,3);
% 
% %lets plot out PV STAs
% axisVal = ceil(sqrt(numPVs));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:numPVs
%     subplot(axisVal,axisVal,i)
%     imagesc(pvStore(:,:,i))
%     colormap('parula')
%     title(num2str(pvSpikeNum(i)))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %plot out MSN STAs
% axisVal = ceil(sqrt(numMSNs));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:numMSNs
%     subplot(axisVal,axisVal,i)
%     imagesc(msnStore(:,:,i))
%     colormap('parula')
%     title(num2str(msnSpikeNum(i)))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %now lets limit by MSNs with greater than 200 spikes
% msns200 = find(msnSpikeNum > 200);
% axisVal = ceil(sqrt(length(msns200)));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:length(msns200)
%     subplot(axisVal,axisVal,i)
%     imagesc(msnStore(:,:,msns200(i)))
%     colormap('parula')
%     title(num2str(msnSpikeNum(msns200(i))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% 
% %% Lets try this while organizing by number of spikes.
% [B I] = sort(pvSpikeNum);
% %lets plot out PV STAs
% axisVal = ceil(sqrt(numPVs));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:numPVs
%     subplot(axisVal,axisVal,i)
%     imagesc(pvStore(:,:,I(i)))
%     colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %plot out MSN STAs
% [B I] = sort(msnSpikeNum);
% axisVal = ceil(sqrt(numMSNs));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:numMSNs
%     subplot(axisVal,axisVal,i)
%     imagesc(msnStore(:,:,I(i)))
%     colormap('parula')
%     title(num2str(msnSpikeNum(I(i))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %now lets limit by MSNs with greater than 200 spikes
% msns200 = find(msnSpikeNum > 200);
% tester = msnSpikeNum(msns200);
% [B I] = sort(tester);
% 
% axisVal = ceil(sqrt(length(msns200)));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:length(msns200)
%     subplot(axisVal,axisVal,i)
%     imagesc(msnStore(:,:,msns200(I(i))))
%     colormap('parula')
%     title(num2str(msnSpikeNum(msns200(I(i)))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %% Lets remove the post-zero points, and order by number of spikes.
% [B I] = sort(pvSpikeNum);
% %lets plot out PV STAs
% axisVal = ceil(sqrt(numPVs));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:numPVs
%     subplot(axisVal,axisVal,i)
%     imagesc(pvStore(:,1:200,I(i)))
%     colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %plot out MSN STAs
% [B I] = sort(msnSpikeNum);
% axisVal = ceil(sqrt(numMSNs));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:numMSNs
%     subplot(axisVal,axisVal,i)
%     imagesc(msnStore(:,1:200,I(i)))
%     colormap('parula')
%     title(num2str(msnSpikeNum(I(i))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end
% 
% %now lets limit by MSNs with greater than 200 spikes
% msns200 = find(msnSpikeNum > 200);
% tester = msnSpikeNum(msns200);
% [B I] = sort(tester);
% 
% axisVal = ceil(sqrt(length(msns200)));
% hFig = figure;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
% set(hFig, 'Position', [10 80 1900 1000])
% for i = 1:length(msns200)
%     subplot(axisVal,axisVal,i)
%     imagesc(msnStore(:,1:200,msns200(I(i))))
%     colormap('parula')
%     title(num2str(msnSpikeNum(msns200(I(i)))))
%     set(gca,'XTick',[], 'YTick', [])
%     set(gca,'YDir','normal')
% end

%% Lets also plot binned spikes and remove the post-zero points, and order by number of spikes.
[B I] = sort(pvSpikeNum);
%lets plot out PV STAs
axisVal = ceil(sqrt(numPVs));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:numPVs
    subplot(axisVal,axisVal*2,2*i-1)
    imagesc(pvStore(:,1:200,I(i)))
    colormap('parula')
    title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    subplot(axisVal,axisVal*2,2*i)
    imagesc(pvTuneStore(:,:,I(i)))
    colormap('parula')
    title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%plot out MSN STAs
[B I] = sort(msnSpikeNum);
axisVal = ceil(sqrt(numMSNs));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:numMSNs
    subplot(axisVal,axisVal*2,2*i-1)
    imagesc(msnStore(:,1:200,I(i)))
    colormap('parula')
    title(num2str(msnSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    subplot(axisVal,axisVal*2,2*i)
    imagesc(msnTuneStore(:,:,I(i)))
    colormap('parula')
    title(num2str(msnSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%now lets limit by MSNs with greater than 200 spikes
msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);

axisVal = ceil(sqrt(length(msns200)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msns200)
    
    subplot(axisVal,axisVal*2,2*i-1)
    imagesc(msnStore(:,1:200,msns200(I(i))))
    colormap('parula')
    title(num2str(msnSpikeNum(msns200(I(i)))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    subplot(axisVal,axisVal*2,2*i)
    imagesc(msnTuneStore(:,:,msns200(I(i))))
    colormap('parula')
    title(num2str(msnSpikeNum(msns200(I(i)))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%% Lets also plot average histogram, binned spikes and remove the post-zero points, and order by number of spikes.
[B I] = sort(pvSpikeNum);
%find zeros and 100 ms from histogram vector
histBinVect = s.(s.DesignationName{1}).HistBinVector;
findZero = find(histBinVect > 0,1,'first');
find100 = find(histBinVect > 0.1,1,'first');
%lets plot out PV STAs
axisVal = ceil(sqrt(numPVs));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:numPVs
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(pvHistStore(:,I(i)),'k')
    plot([findZero findZero],[min(pvHistStore(:,I(i))) max(pvHistStore(:,I(i)))],'b')
    plot([find100 find100],[min(pvHistStore(:,I(i))) max(pvHistStore(:,I(i)))],'b')
    xlim([80 180])
    ylim([min(pvHistStore(:,I(i))) max(pvHistStore(:,I(i)))])
    set(gca,'XTick',[], 'YTick', [])
    title(num2str(pvSpikeNum(I(i))))
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(pvStore(:,1:200,I(i)))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(pvTuneStore(:,:,I(i)))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%plot out MSN STAs
[B I] = sort(msnSpikeNum);
axisVal = ceil(sqrt(numMSNs));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:numMSNs
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(msnHistStore(:,I(i)),'k')
    plot([findZero findZero],[min(msnHistStore(:,I(i))) max(msnHistStore(:,I(i)))],'b')
    plot([find100 find100],[min(msnHistStore(:,I(i))) max(msnHistStore(:,I(i)))],'b')
    xlim([80 180])
    ylim([min(msnHistStore(:,I(i))) max(msnHistStore(:,I(i)))])
    set(gca,'XTick',[], 'YTick', [])
    title(num2str(msnSpikeNum(I(i))))
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,I(i)))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(msnTuneStore(:,:,I(i)))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%now lets limit by MSNs with greater than 200 spikes
msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);

axisVal = ceil(sqrt(length(msns200)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msns200)
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(msnHistStore(:,msns200(I(i))),'k')
    plot([findZero findZero],[min(msnHistStore(:,msns200(I(i)))) max(msnHistStore(:,msns200(I(i))))],'b')
    plot([find100 find100],[min(msnHistStore(:,msns200(I(i)))) max(msnHistStore(:,msns200(I(i))))],'b')
    xlim([80 180])
    ylim([min(msnHistStore(:,msns200(I(i)))) max(msnHistStore(:,msns200(I(i))))])
    set(gca,'XTick',[], 'YTick', [])
%     title(num2str(msnSpikeNum(msns200(I(i)))))
    title(i)
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,msns200(I(i))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(i)))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(msnTuneStore(:,:,msns200(I(i))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(i)))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%% Lets plot out hand selected listings:
[B I] = sort(pvSpikeNum);
pvSel = [1,7,9,12,13,14,17,18,19,21,23,24,25,27,28,30];
%find zeros and 100 ms from histogram vector
histBinVect = s.(s.DesignationName{1}).HistBinVector;
findZero = find(histBinVect > 0,1,'first');
find100 = find(histBinVect > 0.1,1,'first');
%lets plot out PV STAs
axisVal = ceil(sqrt(length(pvSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(pvSel)
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(pvHistStore(:,I(pvSel(i))),'k')
    plot([findZero findZero],[min(pvHistStore(:,I(pvSel(i)))) max(pvHistStore(:,I(pvSel(i))))],'b')
    plot([find100 find100],[min(pvHistStore(:,I(pvSel(i)))) max(pvHistStore(:,I(pvSel(i))))],'b')
    xlim([80 180])
    ylim([min(pvHistStore(:,I(pvSel(i)))) max(pvHistStore(:,I(pvSel(i))))])
    set(gca,'XTick',[], 'YTick', [])
    title(num2str(pvSpikeNum(I(pvSel(i)))))
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(pvStore(:,1:200,I(pvSel(i))))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(pvTuneStore(:,:,I(pvSel(i))))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

hold off
spikeGraphName = 'PVselectedHist-STA-Curve';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets limit by MSNs with greater than 200 spikes
msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);
msnSel = [5,6,7,8,9,10,18,19,22,25,26,29,31,34,35,36,37,39,40,42,43,44,45,48,49,51,54,55,57,58,62,63,67,68,69,70,72,73,74,76,77,78,79,80,83,84,89,90,91,92,93,94,95,96,97,98,102,103,105,106,107,110,111,112,114,115,116,119,120,121];

axisVal = ceil(sqrt(length(msnSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msnSel)
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(msnHistStore(:,msns200(I(msnSel(i)))),'k')
    plot([findZero findZero],[min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))],'b')
    plot([find100 find100],[min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))],'b')
    xlim([80 180])
    ylim([min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))])
    set(gca,'XTick',[], 'YTick', [])
%     title(num2str(msnSpikeNum(msns200(I(msnSel(i))))))
    title(i)
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(msnTuneStore(:,:,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

%lets try splitting in half. 
msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);
msnSel = [5,6,7,8,9,10,18,19,22,25,26,29,31,34,35,36,37,39,40,42,43,44,45,48,49,51,54,55,57,58,62,63,67,68,69];

axisVal = ceil(sqrt(length(msnSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msnSel)
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(msnHistStore(:,msns200(I(msnSel(i)))),'k')
    plot([findZero findZero],[min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))],'b')
    plot([find100 find100],[min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))],'b')
    xlim([80 180])
    ylim([min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))])
    set(gca,'XTick',[], 'YTick', [])
    title(num2str(msnSpikeNum(msns200(I(msnSel(i))))))
%     title(i)
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(msnTuneStore(:,:,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end


hold off
spikeGraphName = 'MSNselectedHist-STA-Curve1';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);
msnSel = [70,72,73,74,76,77,78,79,80,83,84,89,90,91,92,93,94,95,96,97,98,102,103,105,106,107,110,111,112,114,115,116,119,120,121];

axisVal = ceil(sqrt(length(msnSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msnSel)
    if i == 1
        subplot(axisVal*2,axisVal,rem(i,axisVal) + (axisVal)*floor((i)/axisVal))
    else
        subplot(axisVal*2,axisVal,rem(i,axisVal) + axisVal*floor(i/axisVal) + (axisVal)*floor((i-1)/axisVal))
    end
    hold on
    plot(msnHistStore(:,msns200(I(msnSel(i)))),'k')
    plot([findZero findZero],[min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))],'b')
    plot([find100 find100],[min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))],'b')
    xlim([80 180])
    ylim([min(msnHistStore(:,msns200(I(msnSel(i))))) max(msnHistStore(:,msns200(I(msnSel(i)))))])
    set(gca,'XTick',[], 'YTick', [])
    title(num2str(msnSpikeNum(msns200(I(msnSel(i))))))
%     title(i)
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2) + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2-1 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
    if i == 1
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)+1 + (axisVal)*2*(floor((i)/axisVal)+1))
    else
        subplot(axisVal*2,axisVal*2,rem(i,axisVal*2)*2 + axisVal*2*floor(i/(axisVal*2)) + (axisVal)*2*(floor((i)/(axisVal*2))+1) + axisVal*2*floor((i-1)/(axisVal)))
    end
%     subplot(axisVal*2,axisVal*2,2*i+axisVal*2)
    imagesc(msnTuneStore(:,:,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end


hold off
spikeGraphName = 'MSNselectedHist-STA-Curve2';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%% now just plot out the STAs.
[B I] = sort(pvSpikeNum);
pvSel = [1,7,9,12,13,14,17,18,19,21,23,24,25,27,28,30];
%find zeros and 100 ms from histogram vector
histBinVect = s.(s.DesignationName{1}).HistBinVector;
findZero = find(histBinVect > 0,1,'first');
find100 = find(histBinVect > 0.1,1,'first');
%lets plot out PV STAs
axisVal = ceil(sqrt(length(pvSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(pvSel) 

    subplot(axisVal,axisVal,i)
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(pvStore(:,1:200,I(pvSel(i))))
    colormap('parula')
%     title(num2str(pvSpikeNum(I(i))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end

hold off
spikeGraphName = 'PVSelectedSTAs';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


%lets try splitting in half. 
msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);
msnSel = [5,6,7,8,9,10,18,19,22,25,26,29,31,34,35,36,37,39,40,42,43,44,45,48,49,51,54,55,57,58,62,63,67,68,69];

axisVal = ceil(sqrt(length(msnSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msnSel)
    subplot(axisVal,axisVal,i)
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
    
end


hold off
spikeGraphName = 'MSNselectedSTAs1';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

msns200 = find(msnSpikeNum > 200);
tester = msnSpikeNum(msns200);
[B I] = sort(tester);
msnSel = [70,72,73,74,76,77,78,79,80,83,84,89,90,91,92,93,94,95,96,97,98,102,103,105,106,107,110,111,112,114,115,116,119,120,121];

axisVal = ceil(sqrt(length(msnSel)));
hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [10 80 1900 1000])
for i = 1:length(msnSel)
    subplot(axisVal,axisVal,i)
%     subplot(axisVal*2,axisVal*2,2*i-1+axisVal*2)
    imagesc(msnStore(:,1:200,msns200(I(msnSel(i)))))
    colormap('parula')
%     title(num2str(pvSpikeNum(msns200(I(msnSel(i))))))
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','normal')
end


hold off
spikeGraphName = 'MSNselectedSTAs2';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
