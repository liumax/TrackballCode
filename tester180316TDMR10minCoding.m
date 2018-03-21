%This is code to try and make STRFs!



fileName = '180314_ML180306A_R17_3000mid1_dmr10min3x';
% fileName = '180315_ML180306C_R17_3218mid1_3x10minDMR';
%% Constants and things you might want to tweak
%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 0; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 0; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.STRFWin = [-0.4 0.05];
s.Parameters.RPVTime = 0.002; %time limit in seconds for consideration as an RPV
s.Parameters.ClusterWindow = [-0.01 0.03]; %window in seconds for displaying RPV info
s.Parameters.histBin = 0.005; %histogram bin size in seconds
s.Parameters.trodesFS = 30000;%trodes sampling rate

%for duplicate elimination
s.Parameters.DownSampFactor = 10; % how much i want to downsample trodes sampling rate. 3 means sampling every third trodes time point
s.Parameters.corrSlide = 0.05; % window in seconds for xcorr
s.Parameters.ThresholdComparison = 0.05; % percentage overlap to trigger xcorr

%for rotary encoder:
s.Parameters.InterpolationStepRotary = 0.01; %interpolation steps in seconds.

%for edr
s.Parameters.EDRdownsamp = 20; %number of samples to downsample by. Smoothing is likely unnecessary
s.Parameters.EDRTimeCol = 1;
s.Parameters.EDRTTLCol = 3;
s.Parameters.EDRPiezoCol = 2;

%for plotting speed vs firing
s.Parameters.SpeedFiringBins = 1; %bins in seconds for firing rate for display with velocity. 

%for dmr!
s.Parameters.dmrTiming = 6;
s.Parameters.expectedDur = 10; %expected time duration. 
sprFile = 'E:\GIT\cleanDMR\extractedDMR10mindsT5dsF5.mat';
dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
s.Parameters.MinFreq = 4000;
s.Parameters.MaxFreq = 64000;
s.Parameters.DMRtimeDilation = 0.9995; %dilation of time. Multiply Trodes time by this to get DMR time.  
s.Parameters.DMRtimeBin = 1/6000;
newWin = round(s.Parameters.STRFWin/s.Parameters.DMRtimeBin);
preDelay = 0.05;
postDelay = 1/6;

%for plotting the laser related stuff
s.Parameters.LaserWindow = [-0.1 0.1]; %plotting window for rasters
s.Parameters.LaserBin = 0.001; %histogram bin size
s.Parameters.LaserLim = 0.015; %maximum lag value for calculation.

%set other things
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
% format short
disp('Parameters Set')
%% sets up file saving stuff
saveName = strcat(fileName,'FullTuningAnalysis','.mat');
fname = saveName;
pname = pwd;

%% Establishes folders and extracts files!
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';

%% Find and ID Matclust Files for Subsequent Analysis. Generates Structured Array for Data Storage
%pull matclust file names
[matclustFiles] = functionFileFinder(subFoldersCell,'matclust','matclust');
%pull lfp file names. this allows detection of the number of trodes. 
try
    [lfpFiles] = functionFileFinder(subFoldersCell,'LFP','LFP');
    s.NumberTrodes = length(lfpFiles);
catch
    [paramFiles] = functionFileFinder(subFoldersCell,'matclust','param');
    s.NumberTrodes = length(paramFiles) - length(matclustFiles);
end
disp('Matclust Files Extracted')
%generate placeholder structure
% s = struct;
%fill structure with correct substructures (units, not clusters/trodes) and
%then extract waveform and spike data.
[s, truncatedNames] = functionMatclustExtraction(s.Parameters.RPVTime,...
    matclustFiles,s,s.Parameters.ClusterWindow,subFoldersCell);

if toggleDuplicateElimination ==1
    if length(s.DesignationName) > 1
        disp('Now Selecting Based on xCORR')
        [s] = functionDuplicateElimination(s,s.Parameters.DownSampFactor,...
            s.Parameters.corrSlide,s.Parameters.ThresholdComparison,s.Parameters.trodesFS,s.Parameters.RPVTime,s.Parameters.ClusterWindow,...
            s.ShankDesignation,s.ShankMap,s.Shanks);
    end
else
    disp('NOT EXECUTING DUPLICATE ELIMINATION')
end
%pull number of units, as well as names and designation array.
numUnits = size(s.DesignationName,2);
desigNames = s.DesignationName;
desigArray = s.DesignationArray;

%% Extract DIO information. Tuning curve should rely on just one DIO output, DIO1.
disp('Beginning DIO Extraction...')
%find DIO folder and D1 file for analysis
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
if length(D1FileName) == 0
    [D1FileName] = functionFileFinder(subFoldersCell,'DIO','Din1');
end
D1FileName = D1FileName{1};

%extracts DIO stuffs! this code is written to extract inputs for d1
[DIOData] = readTrodesExtractedDataFile(D1FileName);

%pulls out DIO up state onsets.
[dioTimes,dioTimeDiff] = functionBasicDIOCheck(DIOData,s.Parameters.trodesFS);
s.SoundData.ToneTimes = dioTimes;
s.SoundData.ToneTimeDiff = dioTimeDiff;



%pull D2 info (laser ID)
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
if length(D2FileName) == 0
    [D2FileName] = functionFileFinder(subFoldersCell,'DIO','Din2');
end
D2FileName = D2FileName{1};
%extracts DIO stuffs! this code is written to extract inputs for d1
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%pulls out DIO up state onsets.
[dio2Times,dio2TimeDiff] = functionBasicDIOCheck(DIO2Data,s.Parameters.trodesFS);

%check for tons of DIO2. if alot this should be an ID tuning.
if length(dio2Times) > 10
    idToggle = 1;
    disp('Laser Pulses Detected: ID Session')
else
    idToggle = 0;
    disp('No Laser Pulses Detected')
end

%now going into DIO1! Expect a certain binding. 
%watching carefully, it seems like it always misses the first spike. 

%first step, is to isolate the big ITIs. This separates out the chunks 

findBig = find(dioTimeDiff > 10);
findBig(2:end+1) = findBig(1:end);
findBig(1) = 0; %insert first TTL here as well. 
findBig(end+1) = length(dioTimes);

expectedTTLs = s.Parameters.expectedDur*s.Parameters.dmrTiming*60-1; %-1 is fudge factor. 


%% now I want to go through and check on my TTLs: how many are there, are
%they the right ITI, etc. 
badInd = 1;
badStore = [];
badWarn = zeros(length(findBig)-1,1);
badTTLs = [];

whileCounter = 1;
whileTrig = 0;
while whileTrig == 0;
    if whileCounter == length(findBig)-1
        whileTrig = 1;
    end
    targetSamps = dioTimeDiff(findBig(whileCounter)+1:findBig(whileCounter + 1)-1);
    %first confirm the correct number of ITIs. 
    if length(targetSamps) == expectedTTLs
        disp(strcat('Correct Number of TTLs for',num2str(whileCounter)))
        %check if first and second ITI is correct
        if targetSamps(1) > 0.04 & targetSamps(1) < 0.06
            disp(strcat('First Sample Correct for',num2str(whileCounter)))
        else
            error(strcat('First Sample Incorrect for',num2str(whileCounter)))
        end
        if targetSamps(2) > 0.06 & targetSamps(2) < 0.07
            disp(strcat('Second Sample Correct for',num2str(whileCounter)))
        else
            error(strcat('Second Sample Incorrect for',num2str(whileCounter)))
        end
        %now check remainder!
        badFind = find(targetSamps(3:end) > 0.170 | targetSamps(3:end) < 0.16);
        if length(badFind)
            badTTLs = badFind + 2 + findBig(whileCounter)+1;%this will identify the TTL in the dioTimes
            disp(strcat(num2str(length(badTTLs)),' BAD TTLS FOUND'))
            badWarn(whileCounter) = 1;
            badStore(badInd:badInd + length(badTTLs) - 1) = badTTLs;
        end
        whileCounter = whileCounter + 1;
    else
        disp(strcat('Incorrect Number of TTLs for',num2str(whileCounter)))
        disp('Deleting Entire Section Of Recording')
        dioTimes(findBig(whileCounter)+1:findBig(i+1)) = [];
        dioTimeDiff = diff(dioTimes);
        
        findBig = find(dioTimeDiff > 10);
        findBig(2:end+1) = findBig(1:end);
        findBig(1) = 0; %insert first TTL here as well. 
        findBig(end+1) = length(dioTimes);
    end
    
end

trueStarts = dioTimes(findBig(1:end-1)+1) - 0.05*s.Parameters.DMRtimeDilation;
trueEnds = dioTimes(findBig(2:end)-1) + 1/6*s.Parameters.DMRtimeDilation;

timeFilters = zeros(1,2);
timeFilters(1,1) = s.TimeFilterRange(1);
timeFilters(1,2) = trueStarts(1)-s.Parameters.STRFWin(1)+ preDelay;
for i = 2:length(trueStarts)
    timeFilters(i,1) = trueEnds(i-1)-s.Parameters.STRFWin(2) - postDelay;
    timeFilters(i,2) = trueStarts(i)-s.Parameters.STRFWin(1) + preDelay;
end
timeFilters(end+1,1) = trueEnds(end)-s.Parameters.STRFWin(2) - postDelay;
timeFilters(end,2) = s.TimeFilterRange(2);
%now include badTTLs as blocked out times. 
startInd = size(timeFilters,1);
for i = 1:length(badTTLs)
    timeFilters(startInd + i,1) = dioTimes(badTTLs(i)-1)-s.Parameters.STRFWin(2);
    timeFilters(startInd + i,2) = dioTimes(badTTLs(i))-s.Parameters.STRFWin(1);
end

%% Extract data from rotary encoder.
[funcOut] = functionRotaryExtraction(s.Parameters.trodesFS,s.Parameters.InterpolationStepRotary,subFoldersCell);
s.RotaryData = funcOut;
disp('Rotary Encoder Data Extracted')


%% Now I need to load the spr file. 

load(sprFile)

%figure out octave labels
octaveNum = log(s.Parameters.MaxFreq / s.Parameters.MinFreq)/log(2);
numSteps = size(stimulus,1);
stepsPerOct = numSteps/octaveNum;


freqs(1) = s.Parameters.MinFreq;
for i = 2:numSteps
    freqs (i) = freqs(i-1)*(2^(octaveNum/(numSteps-1)));
end


%% Now lets start analysis!
%basic idea is to go through unit by unit and pull out all spikes. Then
%weed out spikes outside of tone presentation times. Then average stimulus
%before that spike from the window specified. Can pull baseline spikes from
%non stimulus periods. 
load(dsTTLFile);


for i = 1:numUnits
    disp(strcat('Working on Unit ',num2str(i)))
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    
    
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    
    
    %weed out bad spike times
    for j = 1:size(timeFilters,1)
        findSpikes = find(spikeTimes >= timeFilters(j,1) & spikeTimes <= timeFilters(j,2));
        spikeTimes(findSpikes) = [];
    end
    averageRate = length(spikeTimes) / (sum(timeFilters(:,2) - timeFilters(:,1)));
    s.(desigNames{i}).AverageRate = averageRate;
    
    
    s.(desigNames{i}).SessionFiring = sessionFiring;
    weedSpikes = [];
    %store weeded spjkes based on DMR number
    for j = 1:length(trueStarts)
        weedSpikes{j} = spikeTimes(spikeTimes < trueEnds(j) & spikeTimes > trueStarts(j));
    end
    
%     spotsPre = round((0.05+0.05+2/30)/s.Parameters.DMRtimeBin);
% %     spotsPre = round((0.05)*s.Parameters.DMRtimeDilation/s.Parameters.DMRtimeBin);
% %     spotsPost = round(1/6*s.Parameters.DMRtimeDilation/s.Parameters.DMRtimeBin);
%     fakeLine = linspace(1,expectedTTLs,length(stimulus)-spotsPre-spotsPost);
    
    %lets pull the stuff based on all spikes
    for j = 1:length(trueStarts)
        disp(strcat('Working on DMR Presentation:',num2str(j)))
%         dmrTimeVector = interp1([1:expectedTTLs-2],([dioTimes(1+(j-1)*(expectedTTLs+1)+2:(j-1)*(expectedTTLs+1)+expectedTTLs)]),fakeLine);
        dmrTimeVector = interp1(highFind,[dioTimes(1+(j-1)*(expectedTTLs+1):(j-1)*(expectedTTLs+1)+expectedTTLs+1)],[highFind(1):highFind(end)]);
%         meanJump = mode(diff(dmrTimeVector));
%         dmrTimeVector(spotsPre+1:end+spotsPre) = dmrTimeVector;
%         fillerBatch = linspace(trueStarts(j),dmrTimeVector(1)-meanJump,spotsPre);
%         dmrTimeVector(1:spotsPre) = fillerBatch;
%         dmrTimeVector(end+1:end+spotsPost) = linspace(dmrTimeVector(end)+meanJump,trueEnds(j),spotsPost);
        targetSpikes = weedSpikes{j};
        indivStore = [];
        for k = 1:length(targetSpikes)
            findTarget = find(dmrTimeVector - targetSpikes(k) > 0,1,'first');
            
            takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
            indivStore(:,:,k) = takeChunk;
        end
        avStore(:,:,j) = mean(indivStore,3);
        bigStore{j} = indivStore;
    end
    bigAvStore(:,:,i) = mean(avStore,3);
%     masterStore{i} = bigStore;
%     masterAvStore{i} = avStore;
    
end



for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(4,6,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    %plots ISI
    subplot(4,6,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
%     %Plot binned response during tone period
%     subplot(4,3,4)
%     imagesc(s.(desigNames{i}).BinTone')
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Mean Binned Response (tone)')
%     %Plot binned response during general period
%     subplot(4,3,7)
%     imagesc(s.(desigNames{i}).BinGen')
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Mean Binned Response (general)')
%     %Plot peak response during general period
%     subplot(4,3,10)
%     imagesc(s.(desigNames{i}).PeakMap')
%     colormap(parula)
%     colorbar
%     set(gca,'XTick',octaveRange(:,2));
%     set(gca,'XTickLabel',octaveRange(:,1));
%     set(gca,'YTick',dbRange(:,2));
%     set(gca,'YTickLabel',dbRange(:,1));
%     title('Peak Response (general)')
    %plot velocity data
    subplot(4,3,6)
    hold on
    plot(s.RotaryData.Velocity(:,1),s.RotaryData.Velocity(:,2)/max(s.RotaryData.Velocity(:,2)),'b')
    plot([s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)],s.(desigNames{i}).SessionFiring/max(s.(desigNames{i}).SessionFiring),'r')
    xlim([s.RotaryData.Velocity(1,1),s.RotaryData.Velocity(end,1)])
    ylim([-0.1,1])
    if toggleROC == 1
        title(strcat('Vel & Firing Rate. AUC:',num2str(s.(desigNames{i}).TrueAUC),'99%Range',num2str(prctile(s.(desigNames{i}).ShuffleAUC,99)),'-',num2str(prctile(s.(desigNames{i}).ShuffleAUC,1))))
    else
        title('Vel & Firing Rate')
    end

%     if idToggle == 0
%         %plot probability of response (tone)
%         subplot(4,3,9)
%         imagesc(s.(desigNames{i}).ProbTone')
%         colormap(parula)
%         colorbar
%         set(gca,'XTick',octaveRange(:,2));
%         set(gca,'XTickLabel',octaveRange(:,1));
%         set(gca,'YTick',dbRange(:,2));
%         set(gca,'YTickLabel',dbRange(:,1));
%         title('Probability of Response (tone)')
%         %plot probability of response (gen)
%         subplot(4,3,12)
%         imagesc(s.(desigNames{i}).ProbGen')
%         colormap(parula)
%         colorbar
%         set(gca,'XTick',octaveRange(:,2));
%         set(gca,'XTickLabel',octaveRange(:,1));
%         set(gca,'YTick',dbRange(:,2));
%         set(gca,'YTickLabel',dbRange(:,1));
%         title('Probability of Response (general)')
%     elseif idToggle == 1
%         %plot laser raster
%         subplot(4,3,9)
%         plot(s.(desigNames{i}).LaserRasters(:,1),s.(desigNames{i}).LaserRasters(:,2),'k.')
%         xlim([s.Parameters.LaserWindow(1),s.Parameters.LaserWindow(2)])
%         title(strcat('Laser Raster, Response %:',num2str(sum(s.(desigNames{i}).LaserResps)/length(dio2Times))))
%         subplot(4,3,12)
%         plot(laserBinVect,s.(desigNames{i}).LaserHist)
%         xlim([s.Parameters.LaserWindow(1),s.Parameters.LaserWindow(2)])
%         title('Laser Histogram')
%     end
% 
%     %plots rasters (chronological)
%     subplot(3,3,2)
%     plot(s.(desigNames{i}).AllRasters(:,1),...
%         s.(desigNames{i}).AllRasters(:,2),'k.','markersize',4)
%     hold on
%     ylim([0 totalTrialNum])
%     xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
%     plot([0 0],[ylim],'b');
%     plot([toneDur toneDur],[ylim],'b');
%     title({fileName;desigNames{i}},'fontweight','bold')
%     set(0, 'DefaulttextInterpreter', 'none')
%     %plots rasters (frequency and amplitude organized)
%     subplot(3,3,5)
%     plot(s.(desigNames{i}).AllRasters(:,1),...
%         s.(desigNames{i}).AllRasters(:,3),'k.','markersize',4)
%     hold on
%     plot([0 0],[ylim],'b');
%     plot([toneDur toneDur],[ylim],'b');
%     rasterFreqLines = zeros(numFreqs,2);
%     if numDBs ==1
%         rasterFreqLines(:,1) = cumsum((matrixTrialNum'));
%     elseif numDBs > 1
%         rasterFreqLines(:,1) = cumsum(sum(matrixTrialNum'));
%     end
% 
%     rasterFreqLines(:,2) = uniqueFreqs;
%     %this generates green lines separating by Frequency
%     tempHold = 1;
%     for k = 1:size(uniqueFreqs,1)
%         plot(s.Parameters.RasterWindow,[tempHold+sum(matrixTrialNum(k,:)) tempHold+sum(matrixTrialNum(k,:))],'g','LineWidth',1)
%         tempHold = tempHold + sum(matrixTrialNum(k,:));
%     end
%     set(gca,'YTick',rasterFreqLines(:,1));
%     set(gca,'YTickLabel',rasterFreqLines(:,2));
%     set(gca,'Ydir','reverse')
%     ylim([0 totalTrialNum])
%     xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
%     title('Descending = increase in amplitude and freq')
%     %plot heatmap organized by frequency
%     subplot(3,3,8)
%     imagesc(s.(desigNames{i}).FrequencyHistograms(:,:))
%     colormap(parula)
%     colorbar
%     set(gca,'YTick',octaveRange(:,2));
%     set(gca,'YTickLabel',octaveRange(:,1));
%     set(gca,'XTick',[1:10:size(histBinVector,2)]);
%     set(gca,'XTickLabel',histBinVector(1:20:end));
%     histBinZero = interp1(histBinVector,1:1:size(histBinVector,2),0);
%     histBinTone = interp1(histBinVector,1:1:size(histBinVector,2),toneDur);
%     line([histBinZero histBinZero],[0 numFreqs],'LineWidth',3,'Color','green')
%     line([histBinZero histBinZero],[0 numFreqs],'LineWidth',2,'Color','black')
%     line([histBinTone histBinTone],[0 numFreqs],'LineWidth',3,'Color','green')
%     line([histBinTone histBinTone],[0 numFreqs],'LineWidth',2,'Color','black')
% %         title('Heatmap by Frequency and Time Max')
%     title('Frequency Arranged Heatmap')
%     % plot histogram.
%     subplot(4,3,3)
%     plot(histBinVector,s.(desigNames{i}).AllHistograms,'k','LineWidth',2)
%     hold on
%     plot(histBinVector,s.(desigNames{i}).AllHistograms - s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
%     plot(histBinVector,s.(desigNames{i}).AllHistograms + s.(desigNames{i}).HistogramStandardDeviation,'b','LineWidth',1)
%     %plot significant values
%     plot(s.(desigNames{i}).AllHistogramSig.Centers(...
%         s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1),...
%         s.(desigNames{i}).AllHistogramSig.Histogram(...
%         s.(desigNames{i}).AllHistogramSig.Histogram(:,3) == 1,1),...
%         'b*')
%     plot(s.(desigNames{i}).AllHistogramSig.Centers(...
%         s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1),...
%         s.(desigNames{i}).AllHistogramSig.Histogram(...
%         s.(desigNames{i}).AllHistogramSig.Histogram(:,4) == 1,1),...
%         'bo')
%     %plot negative values for first tuning
%     plot(s.(desigNames{i}).AllHistogramSig.Centers(...
%         s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1),...
%         s.(desigNames{i}).AllHistogramSig.Histogram(...
%         s.(desigNames{i}).AllHistogramSig.Histogram(:,6) == 1,1),...
%         'k*')
%     plot([0 0],[ylim],'b');
%     plot([toneDur toneDur],[ylim],'b');
%     xlim([s.Parameters.RasterWindow(1) s.Parameters.RasterWindow(2)])
%     title('Histogram')
    subplot(3,3,2)
    imagesc(bigAvStore(:,:,i))
    colormap('parula')
    colorbar

    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-dpdf','-r0')

end







