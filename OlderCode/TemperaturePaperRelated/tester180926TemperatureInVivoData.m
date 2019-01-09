data{1} = data1MW;
data{2} = data3MW;
data{3} = data5MW;
data{4} = data10MW;
data{5} = data15MW;
numStims = 30;
rasterStore = [];
for i = 1:length(data)
    voltDiff = diff(data{i}(:,3));
    bigDiffFind = find(voltDiff > 1);
    bigNegDiff = find(voltDiff < -1);
    if length(bigDiffFind) ~= numStims
        disp('Incorrect number of positive TTLs detected')
        diffDiff = diff(bigDiffFind);
        smallFind = find(diffDiff < 10);
        bigDiffFind(smallFind+1) = [];
    end
    if length(bigNegDiff) ~= numStims
        disp('Incorrect number of negative TTLs detected')
        diffDiff = diff(bigNegDiff);
        smallFind = find(diffDiff < 10);
        bigNegDiff(smallFind+1) = [];
    end
    %now we need to select for the thirty stims. 
    avSub(i) = mean(bigNegDiff - bigDiffFind);
%     meanScaleVal = meanTimeDiff/5;
    
%     meanTimeDiff = 5/avSub(i);
    meanTimeDiff = 0.0031;
    tenSecWindow = 15/meanTimeDiff;
    oneSecWindow = 5/meanTimeDiff;
    raster = [];
    for j = 1:length(bigDiffFind)
        raster(:,j) = data{i}(bigDiffFind(j)-oneSecWindow:bigDiffFind(j)+tenSecWindow,2);
    end
    rasterStore{i} = raster;
    meanResp(:,i) = mean(raster');
    steResp(:,i) = std(raster')/sqrt(30);
    baseMean(i) = mean(meanResp(1614-oneSecWindow:1614,i));


end
%now we need to generate a correct time vector. It is clear that the timing
%is off a little bit, so lets correct for this. 
% meanScaleVal = mean(avSub)*meanTimeDiff/5;
% for i = 1:oneSecWindow
%     respVect(i) = 0 - meanTimeDiff*meanScaleVal*(oneSecWindow-i);
% end
% for i = 1:tenSecWindow
%     respVect(oneSecWindow+i) = 0 - meanTimeDiff*meanScaleVal*(oneSecWindow-i);
% end
respVect = [-5:meanTimeDiff:15];

%sketchy fix
meanResp(1615,:) = mean([meanResp(1614,:);meanResp(1616,:)]);
steResp(1615,:) = mean([steResp(1614,:);steResp(1616,:)]);

figure
hold on
for i = 1:length(data)
    subVals(:,i) = meanResp(:,i)-baseMean(i);
    plot(respVect,smooth(meanResp(:,i)-baseMean(i),31),'Color',[i/length(data) 0 0])
end
xlabel('Time (s)')
xlim([-5 15])
ylabel('Change in Temperature (C)')
title('In Vivo Temperature Change from 532 nm light at 1, 3, 5, 10, 15 mW')
spikeGraphName = strcat('MeasuredChange');
savefig(hFig,spikeGraphName);
% 
% %save as PDF with correct name
% set(hFig,'Units','Inches');
% pos = get(hFig,'Position');
% set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig,spikeGraphName,'-dpdf','-r0')

% plot(meanResp)

% 

figure
hold on
for i = 1:5
plot(smooth(subVals(:,i),31)/max(smooth(subVals(:,i),31)),'Color',[i/5 0 0])
end



% For plotting against temperature model
hFig = figure;
set(hFig, 'Position', [10 80 800 800])
hold on
for i = 1:length(data)
    subVals(:,i) = meanResp(:,i)-baseMean(i);
    plot(respVect,smooth(subVals(:,i),31) + i,'Color',[i/length(data) 0 0])
    plot(respVect,smooth(subVals(:,i),31) + i - smooth(steResp(:,i),31),'Color',[i/length(data) 0 0])
    plot(respVect,smooth(subVals(:,i),31) + i + smooth(steResp(:,i),31),'Color',[i/length(data) 0 0])
end
% plot(t-1,timeStore)
xlabel('Time (s)')
xlim([-5 15])
ylabel('Change in Temperature (C)')
title('In Vivo Temperature Change from 532 nm light at 1, 3, 5, 10, 15 mW')
set(gca,'TickDir','out');

spikeGraphName = strcat('plotMeasuredWithSTE');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')






