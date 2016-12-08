%This code is meant as follow-up plotting for pairing. 

%Things that need to be inputted: fileName

%Variables:
totalTimeBins = 1000;

%% Find analysis file, open in matlab
%Find analysis file
analysisName = strcat(fileName,'PairingAnalysis.mat');

%open analysis file:
s = open(analysisName);
s = s.s;

%% Select only tuned units
%open decision-tuning data
tuningDecision = find(s.decisionTuning==1);
numUnits = length(tuningDecision);

%pull tuned names
tunedNames = s.DesignationName(tuningDecision);

%generate new data array and take only tuned cells into that array
newS = struct;
newS.DesignationName = tunedNames;
for i = 1:numUnits
    newS.(tunedNames{i}) = s.(tunedNames{i});
end

newS.TTLs = s.TTLs;
newS.TimePeriods = s.TimePeriods;
newS.SoundData = s.SoundData;
newS.Parameters = s.Parameters;

%close s. Replace s with newS (for consistency of code)
clear s
s = newS;
clear newS

soundNames = fields(s.SoundData);
histBinVector = [(s.Parameters.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) + ...
    s.Parameters.histBin/2:s.Parameters.histBin:(s.Parameters.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)-s.Parameters.histBin/2];
rasterWindow = [(s.Parameters.rasterWindow(1)*s.SoundData.(soundNames{1}).ToneDur) (s.Parameters.rasterWindow(2)*s.SoundData.(soundNames{1}).ToneDur)];


timeNames = fields(s.TimePeriods);
keyPoints = zeros(length(timeNames),2);
for j=1:length(timeNames)
    keyPoints(j,1) = s.TimePeriods.(timeNames{j})(1);
    keyPoints(j,2) = s.TimePeriods.(timeNames{j})(2);
end
timePoints = unique(keyPoints);
totalTime = timePoints(end) - timePoints(1);
timeBin = totalTime/totalTimeBins;
timeVector = [timePoints(1)+timeBin/2:timeBin:timePoints(end)-timeBin/2];

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.01 0.01]);

%plot out data that I care about!
for i = 1:numUnits
    hFig = figure
    set(hFig, 'Position', [5 5 1280 1000])
    spikeGraphName = strcat(fileName,tunedNames{i},'SpikeAnalysis');
    %plots average waveform
    subplot(4,8,1)
    hold on
    plot(s.(tunedNames{i}).AverageWaveForms(:,2),'LineWidth',2)
    plot(s.(tunedNames{i}).AverageWaveForms(:,1),'r','LineWidth',1)
    plot(s.(tunedNames{i}).AverageWaveForms(:,3),'r','LineWidth',1)
    title({fileName;tunedNames{i};strcat('AverageFiringRate:',num2str(s.(tunedNames{i}).BaselineFiringRate))});
    set(0, 'DefaulttextInterpreter', 'none')
    %plots ISI
    subplot(4,8,2)
    hist(s.(tunedNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(tunedNames{i}).ISIGraph,1000));
    line([s.Parameters.rpvTime s.Parameters.rpvTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.clusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(tunedNames{i}).RPVPercent));...
        strcat(num2str(s.(tunedNames{i}).RPVNumber),'/',num2str(s.(tunedNames{i}).TotalSpikeNumber))});
    %plot out firing rate over the course of the entire recording
    subplot(4,4,5)
    counts = hist(s.(tunedNames{i}).SpikeTimes,timeVector);
    
    hold on
    for k = 1:length(keyPoints)
        rectangle('Position',[keyPoints(k,1) 0 keyPoints(k,2)-keyPoints(k,1) max(counts)],'FaceColor','c','EdgeColor','r')
    end
    plot(timeVector,counts);
    xlim([timePoints(1) timePoints(end)])
    title('Firing Rate over Entire Protocol')
    %next, plot out firing in response to laser
    
    
    %plot out histograms of individual tones. 
    %first, find the number of tones. This assumes no amplitude changes
    uniqueFreqs = s.SoundData.(soundNames{1}).UniqueFreqs;
    numFreqs = length(uniqueFreqs);
    for j = 1:numFreqs
        %Plot for first session
        subplot(numFreqs,4,2+(j-1)*4)
        
        %plot first tuning, before pairing
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{1},'Analysis')).FreqDBHistograms(j,:,:)),'k','LineWidth',2)
        hold on
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{1},'Analysis')).FreqDBHistograms(j,:,:))-...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{1},'Analysis')).FreqDBHistogramErrors(j,:,:)),'k')
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{1},'Analysis')).FreqDBHistograms(j,:,:))+...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{1},'Analysis')).FreqDBHistogramErrors(j,:,:)),'k')
        plot([0 0],[ylim],'b');
        plot([s.SoundData.(soundNames{1}).ToneDur ...
            s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
        xlim([rasterWindow(1) rasterWindow(2)])
        
        %plot second tuning after pairing
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{2},'Analysis')).FreqDBHistograms(j,:,:)),'r','LineWidth',2)
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{2},'Analysis')).FreqDBHistograms(j,:,:))-...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{2},'Analysis')).FreqDBHistogramErrors(j,:,:)),'r')
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{2},'Analysis')).FreqDBHistograms(j,:,:))+...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{2},'Analysis')).FreqDBHistogramErrors(j,:,:)),'r')
        plot([0 0],[ylim],'b');
        plot([s.SoundData.(soundNames{1}).ToneDur ...
            s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
        xlim([rasterWindow(1) rasterWindow(2)])
        if uniqueFreqs(j) == s.SoundData.Pairing1.Frequency
            title(strcat('Histograms for',num2str(uniqueFreqs(j)),'Hz PAIRED'))
        else
            title(strcat('Histograms for',num2str(uniqueFreqs(j)),'Hz'))
        end
        
        
        
        %plot for second session
        subplot(numFreqs,4,3+(j-1)*4)
        
        %plot first tuning, before pairing
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{4},'Analysis')).FreqDBHistograms(j,:,:)),'k','LineWidth',2)
        hold on
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{4},'Analysis')).FreqDBHistograms(j,:,:))-...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{4},'Analysis')).FreqDBHistogramErrors(j,:,:)),'k')
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{4},'Analysis')).FreqDBHistograms(j,:,:))+...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{4},'Analysis')).FreqDBHistogramErrors(j,:,:)),'k')
        plot([0 0],[ylim],'b');
        plot([s.SoundData.(soundNames{1}).ToneDur ...
            s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
        xlim([rasterWindow(1) rasterWindow(2)])
        
        %plot second tuning after pairing
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{5},'Analysis')).FreqDBHistograms(j,:,:)),'r','LineWidth',2)
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{5},'Analysis')).FreqDBHistograms(j,:,:))-...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{5},'Analysis')).FreqDBHistogramErrors(j,:,:)),'r')
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{5},'Analysis')).FreqDBHistograms(j,:,:))+...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{5},'Analysis')).FreqDBHistogramErrors(j,:,:)),'r')
        plot([0 0],[ylim],'b');
        plot([s.SoundData.(soundNames{1}).ToneDur ...
            s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
        xlim([rasterWindow(1) rasterWindow(2)])
        if uniqueFreqs(j) == s.SoundData.Pairing1.Frequency
            title(strcat('Histograms for',num2str(uniqueFreqs(j)),'Hz PAIRED'))
        else
            title(strcat('Histograms for',num2str(uniqueFreqs(j)),'Hz'))
        end
        
        
        %plot for third session
        subplot(numFreqs,4,4+(j-1)*4)
        
        %plot first tuning, before pairing
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{7},'Analysis')).FreqDBHistograms(j,:,:)),'k','LineWidth',2)
        hold on
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{7},'Analysis')).FreqDBHistograms(j,:,:))-...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{7},'Analysis')).FreqDBHistogramErrors(j,:,:)),'k')
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{7},'Analysis')).FreqDBHistograms(j,:,:))+...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{7},'Analysis')).FreqDBHistogramErrors(j,:,:)),'k')
        plot([0 0],[ylim],'b');
        plot([s.SoundData.(soundNames{1}).ToneDur ...
            s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
        xlim([rasterWindow(1) rasterWindow(2)])
        
        %plot second tuning after pairing
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{8},'Analysis')).FreqDBHistograms(j,:,:)),'r','LineWidth',2)
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{8},'Analysis')).FreqDBHistograms(j,:,:))-...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{8},'Analysis')).FreqDBHistogramErrors(j,:,:)),'r')
        plot(histBinVector,squeeze(s.(tunedNames{i}).(strcat(soundNames{8},'Analysis')).FreqDBHistograms(j,:,:))+...
            squeeze(s.(tunedNames{i}).(strcat(soundNames{8},'Analysis')).FreqDBHistogramErrors(j,:,:)),'r')
        plot([0 0],[ylim],'b');
        plot([s.SoundData.(soundNames{1}).ToneDur ...
            s.SoundData.(soundNames{1}).ToneDur],[ylim],'b');
        xlim([rasterWindow(1) rasterWindow(2)])
        if uniqueFreqs(j) == s.SoundData.Pairing1.Frequency
            title(strcat('Histograms for',num2str(uniqueFreqs(j)),'Hz PAIRED'))
        else
            title(strcat('Histograms for',num2str(uniqueFreqs(j)),'Hz'))
        end
        
    end
    
end







