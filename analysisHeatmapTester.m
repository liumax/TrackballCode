
tester = matclustStruct.matclust_param_nt1.Rasters{1,4};

ampHold = tester(:,4);
[B,I] = sort(ampHold);
ampSort = tester(I,:);
freqHold =ampSort(:,3);
[B,I] = sort(freqHold);
ampFreqSort = ampSort(I,:);

backup = ampFreqSort;
%This code will replace the actual trial number with a flat number for all
%of the same amplitude and dB. For example, all trials with 4kHz at 50dB
%will be labeled with a single number. This collapses individual trials.
for i = 1:9 %this is for frequencies
    for j = 1:6 %for amplitude
        findSeg = find(ampFreqSort(:,3) == uniqueFreqs(i) & ampFreqSort(:,4) == uniqueDBs(j));
        ampFreqSort(findSeg,1) = (6*(i-1))+j;
    end
end

%plotting code that includes horizontal lines at intervals of 6 (to
%separate frequencies)
figure
plot(ampFreqSort(:,2),ampFreqSort(:,1),'k.')
hold on
plot([-0.5 0.5],[6 6])
plot([-0.5 0.5],[12 12])
plot([-0.5 0.5],[18 18])
plot([-0.5 0.5],[24 24])
plot([-0.5 0.5],[30 30])
plot([-0.5 0.5],[36 36])
plot([-0.5 0.5],[42 42])
plot([-0.5 0.5],[48 48])
plot([-0.5 0.5],[54 54])
plot([-0.5 0.5],[60 60])



%This code is to generate histograms of the responses. freqAmpHist will
%hold the individual histograms for each frequency/amplitude. heatMapTester
%will hold all frequencies vs all amplitudes. 
freqAmpHist = cell(9,6);
heatMapTester = zeros(54,40);

for i = 1:9 %this is for frequencies
    for j = 1:6 %for amplitude
        findSeg = find(ampFreqSort(:,3) == uniqueFreqs(i) & ampFreqSort(:,4) == uniqueDBs(j));
        tempHolder = ampFreqSort(findSeg,2);
        [counts centers] = hist(tempHolder,histBinVector);
        countSize = size(counts);
        centerSize = size(centers);
        if countSize(1)>countSize(2)
            counts = counts';
        end
        if centerSize(1)>centerSize(2)
            centers = centers';
        end
        freqAmpHist{i,j} = [counts'*(1/histBin)/20,centers'];
        heatMapTester((6*(i-1))+j,:) = [counts*(1/histBin)/20];
    end
end

figure
hold on
for i = 1:9
    for j = 1:6
        plot(freqAmpHist{i,j}(:,2),freqAmpHist{i,j}(:,1)/20 + (6*(i-1))+j,'r')
    end
end

plot([-0.5 0.5],[6 6])
plot([-0.5 0.5],[12 12])
plot([-0.5 0.5],[18 18])
plot([-0.5 0.5],[24 24])
plot([-0.5 0.5],[30 30])
plot([-0.5 0.5],[36 36])
plot([-0.5 0.5],[42 42])
plot([-0.5 0.5],[48 48])
plot([-0.5 0.5],[54 54])
plot([-0.5 0.5],[60 60])

figure
imagesc(heatMapTester)
colormap hot
hold on
plot([-0.5 0.5],[6 6])
plot([-0.5 0.5],[12 12])
plot([-0.5 0.5],[18 18])
plot([-0.5 0.5],[24 24])
plot([-0.5 0.5],[30 30])
plot([-0.5 0.5],[36 36])
plot([-0.5 0.5],[42 42])
plot([-0.5 0.5],[48 48])
plot([-0.5 0.5],[54 54])
plot([-0.5 0.5],[60 60])

%this finds all the positions of the unique trial numbers, in the case
%where trials have been sorted by amplitude and frequency.
C = unique(backup(:,1),'stable');

CSort = sort(C);

backupRep = backup(:,1);

for i =1:length(C)
    backupRep(backupRep == CSort(i)) = find(C==CSort(i));
end

backup(:,1) = backupRep;

figure
plot(backup(:,2),backup(:,1),'k.')

freqChange = diff(backup(:,3));
ampChange = diff(backup(:,4));

%This code below is meant to graph rasters in order of frequency and
%amplitude. 

figure
plot(backup(:,2),backup(:,1),'k.')
hold on
plot([-0.5 0.5],[backup(136,1) backup(136,1)])
plot([-0.5 0.5],[backup(294,1) backup(294,1)])
plot([-0.5 0.5],[backup(515,1) backup(515,1)])
plot([-0.5 0.5],[backup(839,1) backup(839,1)])
plot([-0.5 0.5],[backup(1573,1) backup(1573,1)])
plot([-0.5 0.5],[backup(1829,1) backup(1829,1)])
plot([-0.5 0.5],[backup(1945,1) backup(1945,1)])
plot([-0.5 0.5],[backup(2072,1) backup(2072,1)])