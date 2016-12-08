
% cd E:\160907PulseSoundRecording\160902_ML160805B_L_1800_1000msLongRecnoLaserPul
% data = matclustStruct.matclust_param_nt4;
% rasterData = data.Rasters{1};
% 
% histHolder = zeros(301,15);
% for i = 1:15
% tempHolder = rasterData(rasterData(:,5) == i,2);
% histHolder(:,i) = smooth(hist(tempHolder,[-0.1:0.001:0.2]),5);
% rasterHolder{i} = rasterData(rasterData(:,5) == i,2);
% end
% 
% 
% figure
% for i = 1:15
% subplot(1,15,i)
% plot(histHolder(:,i))
% end
% 
% histMax = max(max(histHolder));
% normHist = histHolder/histMax*4;
% for i = 1:15
%     normHist(:,i) = normHist(:,i) + 1*i;
% end
% 
% figure
% plot([-0.1:0.001:0.2],normHist)
% xlim([-0.1 0.2])
% ylim([0 17])
% 

cd E:\160915pulseLaser\160909_ML160718B_R_2200_8kHz80dBPulseWithTermPulsedLas
% data = matclustStruct.nt28;
data = matclustStruct.nt18;
rasterData = data.Rasters{1};
rasterData(:,6) = mod(rasterData(:,6),2);

rasterHolder = cell(10,2);
histHolder = zeros(401,10,2);

for j = 1:2
    for i = 1:10
        tempHolder = rasterData(rasterData(:,5) == i & rasterData(:,6) == j-1,2);
        tempHolder = tempHolder(tempHolder < 0.3 & tempHolder > -0.1);
        histHolder(:,i,j) = smooth(hist(tempHolder,[-0.1:0.001:0.3]),11);
        rasterHolder{i,j} = tempHolder;
    end
end

histMax = squeeze(max(max(histHolder)));

normHist = histHolder;
normHist(:,:,1) = histHolder(:,:,1)/histMax(1);
normHist(:,:,2) = histHolder(:,:,2)/histMax(2);

for i = 1:10
    normHist(:,i,:) = normHist(:,i,:) + 1*i;
end

figure
plot([-0.1:0.001:0.3],normHist(:,:,1),'k')
hold on
plot([-0.1:0.001:0.3],normHist(:,:,2),'r')
ylim([1 11])
xlim([-0.1 0.3])



