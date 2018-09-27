data{1} = data1MW;
data{2} = data3MW;
data{3} = data5MW;
data{4} = data10MW;
data{5} = data15MW;

rasterStore = [];
for i = 1:length(data)
    voltDiff = diff(data{i}(:,3));
    bigDiffFind = find(voltDiff > 1);
    
    meanTimeDiff = mean(diff(data{i}(:,1)));
    tenSecWindow = 10/meanTimeDiff;
    oneSecWindow = 1/meanTimeDiff;
    raster = [];
    for j = 1:length(bigDiffFind)
        raster(:,j) = data{i}(bigDiffFind(j)-oneSecWindow:bigDiffFind(j)+tenSecWindow,2);
    end
    rasterStore{i} = raster;
    meanResp(:,i) = mean(raster');
    baseMean(i) = mean(meanResp(1:oneSecWindow,i));


end
respVect = [-1:meanTimeDiff:10];
figure
hold on
for i = 1:length(data)
    subVals(:,i) = meanResp(:,i)-baseMean(i);
    plot(respVect,smooth(meanResp(:,i)-baseMean(i),91),'Color',[i/length(data) 0 0])
end
xlabel('Time (s)')
xlim([-1 10])
ylabel('Change in Temperature (C)')
title('In Vivo Temperature Change from 473 nm light at 1, 3, 5, 10, 15 mW')
% plot(meanResp)

% 

figure
hold on
for i = 1:5
plot(smooth(subVals(:,i),101)/max(smooth(subVals(:,i),101)),'Color',[i/5 0 0])
end