%this code is to link LFP with timing of tone inputs. 
viewWindow = [-0.5,0.5];

lfp = readTrodesExtractedDataFile('160405_ML160218B_R17_2372_fullTune.LFP_nt1ch1.dat');
dio = readTrodesExtractedDataFile('160405_ML160218B_R17_2372_fullTune.dio_D1.dat');
%counts number of LFP samples
lfpSamples = size(lfp.fields.data,1);
%makes LFP time points based on # of samples and decimation
lfpTimes = (((0:1:lfpSamples-1)*lfp.decimation)+lfp.first_timestamp)';
lfpSignals = lfp.fields.data;
%pulls dio data for both state and time
dioHolder = zeros(size(dio.fields(1).data,1),2);
dioHolder(:,1) = dio.fields(1).data;
dioHolder(:,2) = dio.fields(2).data;
%eliminates non- upswing datapoints
dioHolder = dioHolder(dioHolder(:,2) == 1);
%calculates number of samples in the viewing window
viewSamples = (viewWindow(2)-viewWindow(1))*lfp.clockrate/lfp.decimation;

%holds all LFP traces
lfpHolder = zeros(size(dioHolder,1),viewSamples);
for i = 1:size(dioHolder,1)
    finder = find(lfpTimes>dioHolder(i)+viewWindow(1)*lfp.clockrate,1);
    lfpHolder(i,:) = lfpSignals(finder:finder+viewSamples-1);
end

tester = mean(lfpHolder);
plot(tester);

master(:,1) = soundData.Frequencies;
master(:,2) = soundData.dBs;

uniqueFreqs = unique(master(:,1));
uniqueDBs = unique(master(:,2));
trials = soundData.ToneRepetitions;

freqAmpHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1));

for i = 1:size(uniqueFreqs,1)
    for j = 1:size(uniqueDBs,1)
        trialFinder = intersect(find(master(:,1) == uniqueFreqs(i)),...
            find(master(:,2) == uniqueDBs(j)));
        freqAmpHolder{i,j} = lfpHolder(trialFinder,:);
    end
end

%this is to plot all mean traces by frequency
figure
hold on
for i =1:13
plot(mean(freqAmpHolder{i,6}),'color',rand(3,1))
end