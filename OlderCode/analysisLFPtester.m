%this code is to link LFP with timing of tone inputs. 
lfpWindow = [-0.5,0.5];

lfp = readTrodesExtractedDataFile('160429_bathTestingLFP.LFP_nt1ch1.dat');
dio = readTrodesExtractedDataFile('timingTester.dio_D1.dat');
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
viewSamples = (lfpWindow(2)-lfpWindow(1))*lfp.clockrate/lfp.decimation;

%holds all LFP traces
lfpHolder = zeros(size(dioHolder,1),viewSamples);
for i = 1:size(dioHolder,1)
    finder = find(lfpTimes>dioHolder(i)+lfpWindow(1)*lfp.clockrate,1);
    lfpHolder(i,:) = lfpSignals(finder:finder+viewSamples-1);
end

tester = mean(lfpHolder);
plot(tester);

master(:,1) = soundData.Frequencies;
% master(:,2) = soundData.dBs;

uniqueFreqs = unique(master(:,1));
% uniqueDBs = unique(master(:,2));
trials = soundData.ToneRepetitions;

% freqAmpHolder = cell(size(uniqueFreqs,1),size(uniqueDBs,1));
freqAmpHolder = cell(size(uniqueFreqs,1),1);

% for i = 1:size(uniqueFreqs,1)
%     for j = 1:size(uniqueDBs,1)
%         trialFinder = intersect(find(master(:,1) == uniqueFreqs(i)),...
%             find(master(:,2) == uniqueDBs(j)));
%         freqAmpHolder{i,j} = lfpHolder(trialFinder,:);
%     end
% end

for i = 1:size(uniqueFreqs,1)
    trialFinder = find(master(:,1) == uniqueFreqs(i));
    freqAmpHolder{i} = lfpHolder(trialFinder,:);
end

%this is to plot all mean traces by frequency
% figure
% hold on
% for i =length(uniqueFreqs)
% plot(mean(freqAmpHolder{i}),'color',rand(3,1))
% end

h = figure;
hold on
for i =1:length(uniqueFreqs)
plot(mean(freqAmpHolder{i}),'color',rand(3,1))
end
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')