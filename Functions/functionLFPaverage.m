function [lfpStruct] = functionLFPaverage(master, lfpWindow, s,homeFolder,fileName, uniqueFreqs, uniqueDBs, numFreqs, numDBs);

%%
%now I will extract LFP information!
disp('Starting LFP Analysis')
%this code picks out the LFP folder and moves to that folder
lfpFinder =dir;
lfpFinder = {lfpFinder.name};
lfpIndex = strfind(lfpFinder,'LFP');
lfpIndex = find(not(cellfun('isempty', lfpIndex)));
lfpFolder = lfpFinder{lfpIndex};
newDir = strcat(pwd,'\',lfpFolder);
cd(newDir)
%this code pulls the file names of the LFP files, and reorders them in
%natural order (1,2,10, not 1,10,2)
lfpFinder = dir;
lfpFinder = {lfpFinder.name};
lfpIndex = strfind(lfpFinder,'LFP');
lfpIndex = find(not(cellfun('isempty', lfpIndex)));
lfpFiles = cell(size(lfpIndex,2),1);
lfpFiles= lfpFinder(lfpIndex);
[cs index] = sort_nat(lfpFiles);
lfpFiles = cs;
numLFPs = size(lfpFiles,2);
%this generates the set of colors I want to use
colorArray = zeros(60,3);
colorArray(1,:) = [0,0,1];
for lfpInd1 = 1:29
    colorArray(lfpInd1+1,:) = [0,(lfpInd1)/29,(29-lfpInd1)/29];
end
for lfpInd1 = 1:30
    colorArray(lfpInd1+30,:) = [(lfpInd1)/30,(30-lfpInd1)/30,0];
end 
%this spaces out my frequencies to tile the color space I have defined
spacer = size(colorArray,1)/size(uniqueFreqs,1);
spacerArray = 1:1:size(uniqueFreqs,1);
spacerArray = round(spacerArray*spacer);
plotColorsFreq = colorArray(spacerArray,:);

%this spaces things out in the same way by DB
spacer = size(colorArray,1)/size(uniqueDBs,1);
spacerArray = 1:1:size(uniqueDBs,1);
spacerArray = round(spacerArray*spacer);
plotColorsDB = colorArray(spacerArray,:);

%next, I will need to find all the times for both target frequencies and
%target DB ratings. This will provide me with the times with which I need
%to align things!
numTrials = size(master(:,1),1);
targetTimesFreqs = zeros(numTrials/numFreqs,numFreqs);
targetTimesDBs = zeros(numTrials/numDBs,numDBs);

for lfpInd1 = 1:numFreqs
    freqFinder = find(master(:,2) == uniqueFreqs(lfpInd1));
    targetTimesFreqs(:,lfpInd1) = master(freqFinder,1);
end

for lfpInd1 = 1:numDBs
    dbFinder = find(master(:,3) == uniqueDBs(lfpInd1));
    targetTimesDBs(:,lfpInd1) = master(dbFinder,1);
end

%the next step is to extract the specific window I'm interested in from the
%LFP data, and then average it for the given metric (frequency or db)

%first i figure out how many samples are in the LFP analysis window:
lfp = readTrodesExtractedDataFile(lfpFiles{1});
viewSamples = round((lfpWindow(2)-lfpWindow(1))*lfp.clockrate/lfp.decimation);

lfpAverageFreqs = zeros(numLFPs,numFreqs,viewSamples);

for lfpInd1 = 1:numLFPs
        %opens LFP file
    lfp = readTrodesExtractedDataFile(lfpFiles{lfpInd1});
    %counts number of LFP samples
    lfpSamples = size(lfp.fields.data,1);
    %makes LFP time points based on # of samples and decimation
    %adjust time to actual time (to match master)
    lfpTimes = (((0:1:lfpSamples-1)*lfp.decimation)+lfp.first_timestamp)'/30000;
    lfpSignals = lfp.fields.data;
    %generates a temporary holder for LFP traces
    tempHolder = zeros(numTrials/numFreqs,viewSamples);
    %goes through all trials of a given frequency, extracts LFP signals
    %and places into tempHolder
    for lfpInd2 = 1:numFreqs
        for lfpInd3 = 1:numTrials/numFreqs
            finder = find(lfpTimes > targetTimesFreqs(lfpInd3,lfpInd2)+lfpWindow(1),1);
            tempHolder(lfpInd3,:) = lfpSignals(finder:finder + viewSamples - 1);
        end
        lfpAverageFreqs(lfpInd1,lfpInd2,:) = mean(tempHolder);
    end
    disp(strcat('Averaged Electrode ',lfpInd1,' By Frequency'))
end
disp('Done with LFPs by Freq')

lfpAverageDBs = zeros(numLFPs,numDBs,viewSamples);

for lfpInd1 = 1:numLFPs
    %opens LFP file
    lfp = readTrodesExtractedDataFile(lfpFiles{lfpInd1});
    %counts number of LFP samples
    lfpSamples = size(lfp.fields.data,1);
    %makes LFP time points based on # of samples and decimation
    %adjust time to actual time (to match master)
    lfpTimes = (((0:1:lfpSamples-1)*lfp.decimation)+lfp.first_timestamp)'/30000;
    lfpSignals = lfp.fields.data;
    %generates a temporary holder for LFP traces
    tempHolder = zeros(numTrials/numDBs,viewSamples);
    %goes through all trials of a given frequency, extracts LFP signals
    %and places into tempHolder
    for lfpInd2 = 1:numDBs
        for lfpInd3 = 1:numTrials/numDBs
            finder = find(lfpTimes > targetTimesDBs(lfpInd3,lfpInd2)+lfpWindow(1),1);
            tempHolder(lfpInd3,:) = lfpSignals(finder:finder + viewSamples - 1);
        end
        lfpAverageDBs(lfpInd1,lfpInd2,:) = mean(tempHolder);
    end
    disp(strcat('Averaged Electrode ',lfpInd1,' By DB'))
end

disp('Done with LFPs by DB')

%calculate means for LFPs, for both frequency and DB. This is to set y
%limits on graphs

minLFPfreq = min(min(min(lfpAverageFreqs)));
maxLFPfreq = max(max(max(lfpAverageFreqs)));

minLFPdb = min(min(min(lfpAverageDBs)));
maxLFPdb = max(max(max(lfpAverageDBs)));

%calculates the zero point based on number of samples. This is to generate
%lines that indicate tone duration
totalLFPWindow = lfpWindow(2)-lfpWindow(1);
lfpZero = abs(lfpWindow(1))*viewSamples/totalLFPWindow;
toneEnd = abs(s.SoundData.ToneDuration)*viewSamples/totalLFPWindow;

%save everything to a dummy structured array
lfpStruct.MeansByFreq = lfpAverageFreqs;
lfpStruct.MeansByDB = lfpAverageDBs;
lfpStruct.AnalysisWindow = lfpWindow;
lfpStruct.ToneStartTime = lfpZero;
lfpStruct.ToneEndTime = toneEnd;

%Plotting!!

%generates cell array of frequency names for use in legend
freqNameHolder = cell(1,numFreqs);
for lfpInd1 =1:numFreqs
    freqNameHolder{lfpInd1} = num2str(uniqueFreqs(lfpInd1));
end

%generates cell array of DB names for use in legend
dbNameHolder = cell(1,numDBs);
for lfpInd1 =1:numDBs
    dbNameHolder{lfpInd1} = num2str(uniqueDBs(lfpInd1));
end

h = figure;
set(h, 'Position', [100 100 1000 1000])
%generates text box with information!
mTextBox = uicontrol('style','text');
descr = {'LFP By Frequency and DB';
    'File:';
    fileName;
    'Window (s):';
    lfpWindow;
    'Tone Duration(s):';
    s.SoundData.ToneDuration};
set(mTextBox,'String',descr);
set(mTextBox,'Units','normalized');
set(mTextBox,'Position',[0.36,0.74,0.3,0.2])

%plots by frequency
for lfpInd1 = 1:numLFPs
    subplot(numLFPs,3,1+(3*(lfpInd1-1)))
    plot([lfpZero lfpZero],[minLFPfreq maxLFPfreq],'k');
    hold on
    plot([lfpZero+toneEnd lfpZero+toneEnd],[minLFPfreq maxLFPfreq],'k')
    for lfpInd2 = 1:numFreqs
        plot(squeeze(lfpAverageFreqs(lfpInd1,lfpInd2,:)),'color',plotColorsFreq(lfpInd2,:));
        freqTicker(lfpInd2) = plot(squeeze(lfpAverageFreqs(lfpInd1,lfpInd2,:)),'color',plotColorsFreq(lfpInd2,:));
    end
    
    xlim([0 viewSamples])
    ylim([minLFPfreq maxLFPfreq])
    set(gca, 'XTick', [], 'YTick', [])
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1 1.4]) % stretch its width and height
    
    subplot(numLFPs,3,3+(3*(lfpInd1-1)))
    plot([lfpZero lfpZero],[minLFPdb maxLFPdb],'k');
    hold on
    plot([lfpZero+toneEnd lfpZero+toneEnd],[minLFPdb maxLFPdb],'k')
    hold on
    for lfpInd2 = 1:numDBs
        plot(squeeze(lfpAverageDBs(lfpInd1,lfpInd2,:)),'color',plotColorsDB(lfpInd2,:));
        dbTicker(lfpInd2) = plot(squeeze(lfpAverageDBs(lfpInd1,lfpInd2,:)),'color',plotColorsDB(lfpInd2,:));
    end
    
    xlim([0 viewSamples])
    ylim([minLFPdb maxLFPdb])
    set(gca, 'XTick', [], 'YTick', [])
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1 1.4]) % stretch its width and height
end

%generates legend for frequencies
hL = legend(freqTicker,freqNameHolder);
set(hL,'Position', [0.36 0.4 0.3 0.2],'Units','normalized');

hL2 = legend(dbTicker,dbNameHolder);
set(hL2,'Position', [0.36 0.15 0.3 0.2],'Units','normalized');

%returns to original directory
cd(homeFolder)

%save as matlab figure with correct name (fileName+LFP)
lfpName = strcat(fileName,'LFPGraph');
savefig(h,lfpName);

%save as PDF with correct name
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,lfpName,'-dpdf','-r0')
