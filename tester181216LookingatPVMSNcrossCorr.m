cd Z:\Max\181002LookingAtDMR\180718_ML180619B_L_AudStr_pen1_3000_3x10minDMRttlFix
msnData = s.nt14cluster3.SpikeTimes;
fsiData = s.nt13cluster1.SpikeTimes;

msnData = s.nt16cluster3.SpikeTimes;
fsiData = s.nt16cluster4.SpikeTimes;



counter = 1;
rasterStore = [];
for i = 1:length(fsiData)
spikeSub = msnData - fsiData(i);
findSpikes = find(spikeSub > -0.01 & spikeSub < 0.01);
foundSpikes = spikeSub(findSpikes);
if foundSpikes
rasterStore(counter:counter+length(foundSpikes)-1,1) = foundSpikes;
rasterStore(counter:counter+length(foundSpikes)-1,2) = i;
counter = counter + length(foundSpikes);
end
end
figure
hist(rasterStore(:,1),[-0.01:0.0005:0.01])