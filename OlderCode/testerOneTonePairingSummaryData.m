%load s file
load('161104_ML160921C_L17_2800_secondOneTonePairingPairingAnalysis.mat')

%pull significant spikes
sig1 = s.SignificantSpikes.spikesTuning1;
sig2 = s.SignificantSpikes.spikesTuning2;

sigs = sig1+sig2;

sigFinder = find(sigs > 0);
%pull names of significant things
desigNames = s.DesignationName;
% desigNames = desigNames(sigFinder);


%I know I want windows 2 and 3
wantWindows = [2:3];

%get the number of frequencies
numFreqs = length(s.SoundData.Tuning1.UniqueFreqs);

%pull data from individual cells

binDataHolder = zeros(length(sigFinder),numFreqs,4);

for i = 1:length(sigFinder)
    binDataHolder(i,:,1) = mean(squeeze(s.(desigNames{sigFinder(i)}).Tuning1Analysis.BinSpikeStats(:,1,1,wantWindows)),2);
    binDataHolder(i,:,2) = mean(squeeze(s.(desigNames{sigFinder(i)}).Tuning2Analysis.BinSpikeStats(:,1,1,wantWindows)),2);
    binDataHolder(i,:,3) = s.(desigNames{sigFinder(i)}).Tuning1Analysis.BinSpikeStats(:,1,1,1);
    binDataHolder(i,:,4) = s.(desigNames{sigFinder(i)}).Tuning2Analysis.BinSpikeStats(:,1,1,1);
end

subBinDataHolder = binDataHolder(:,:,[1:2]) - binDataHolder(:,:,[3:4]);

binMax = max(max(max(binDataHolder)));

figure
hold on
plot(binDataHolder(:,1,1),binDataHolder(:,1,2),'k.')
plot(binDataHolder(:,2,1),binDataHolder(:,2,2),'b.')
plot(binDataHolder(:,3,1),binDataHolder(:,3,2),'g.')
plot(binDataHolder(:,4,1),binDataHolder(:,4,2),'r.')
plot([0 binMax],[0 binMax],'k')
xlim([0 binMax])
ylim([0 binMax])


figure
plot(mean(binDataHolder(:,:,1),2),mean(binDataHolder(:,:,2),2),'k.')
hold on
plot([0 binMax],[0 binMax],'k')


ratio = mean(binDataHolder(:,:,1),2)./mean(binDataHolder(:,:,2),2);


