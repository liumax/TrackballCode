%This is code to try and make STRFs!


fileName = '180718_ML180619B_L_AudStr_pen1_3000_3x10minDMRttlFix';
% fileName = '180622_ML180515C_R17_pen3_2640_3x10minDMRTTLADJUST';
% fileName = '180315_ML180306C_R17_3218mid1_3x10minDMR';
%% Constants and things you might want to tweak
%TOGGLES FOR ENABLING/DISABLING FEATURES
s.Parameters.toggleRPV = 0; %1 means you use RPVs to eliminate units. 0 means not using RPVs
toggleTuneSelect = 0; %1 means you want to select tuning manually, 0 means no selection.
toggleDuplicateElimination = 0; %1 means you want to eliminate duplicates.
toggleROC = 0; %toggle for tuning on/off ROC analysis

%PARAMETERS FOR BASIC ARRANGEMENT OF DATA
s.Parameters.STRFWin = [-0.1 0.01];
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
% dsTTLFile = 'E:\GIT\cleanDMR\DMR10mindsT5ttlPos.mat';
ttlInfo = 'Z:\Max\dmrOutputFiles\dmr-400flo-64000fhi-4SM-40TM-40db-192000khz-6DF-10min40dBTTLADJUSTdmrStimTimeAndTTLTimes.mat';
s.Parameters.MinFreq = 4000;
s.Parameters.MaxFreq = 64000;
% s.Parameters.DMRtimeDilation = 0.9995; %dilation of time. Multiply Trodes time by this to get DMR time.  
% s.Parameters.DMRtimeBin = 1/6000;
% newWin = round(s.Parameters.STRFWin/s.Parameters.DMRtimeBin);
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

%% QC of TTLs

%first step, is to isolate the big ITIs. This separates out the chunks 

findBig = find(dioTimeDiff > 10);
ttlChunks = zeros(length(findBig)+1,2);
for i = 1:length(ttlChunks)
    if i ==1
        ttlChunks(1,1) = 1;
        ttlChunks(1,2) = findBig(1);
    elseif i <= length(findBig)
        ttlChunks(i,1) = findBig(i-1) + 1;
        ttlChunks(i,2) = findBig(i);
    else
        ttlChunks(i,1) = findBig(i-1) + 1;
        ttlChunks(i,2) = length(dioTimes);
    end
end

%now lets determine if the length of any one segment is incorrect
blockLengths = ttlChunks(:,2) - ttlChunks(:,1)+1;
%load the data file.
load(ttlInfo)
lengthTrue = length(ttlOnsetTime);
for i = 1:length(blockLengths)
    if blockLengths(i) == lengthTrue
        disp(strcat('Block-',num2str(i),'-has correct number of TTLs'))
    else
        error(strcat('Block-',num2str(i),'-has incorrect number of TTLs'))
    end
end
%now lets go through more carefully in each block to look at the itis of
%these TTL pulses
devCounter = 1; %this counts number of deviating TTLs
devStore = [];
for i = 1:length(blockLengths)
    %display block number
    disp(strcat('Crawling through block-',num2str(i)))
    %now pull up appropriate ITIs
    blockITIs = diff(dioTimes(ttlChunks(i,1):ttlChunks(i,2)));
    %compare against true values with subtraction
    itiMod = ttlDiffTime - blockITIs;
    %find all remainders that are too big
    bigMod = find(abs(itiMod) > 0.001);
    if bigMod
        disp(strcat(num2str(length(bigMod)),'-deviations found, marking...'))
        for j = 1:length(bigMod)
            devStore(devCounter,1) = i; %store the chunk
            devStore(devCounter,2) = bigMod(i) - 1; %stores the TTLs before and after. 
            devStore(devCounter,3) = bigMod(i) +1;
            devStore(devCounter,4) = ttlChunks(i,1) - 1 + bigMod(i) - 1; %stores the TTLs before and after. 
            devStore(devCounter,5) = ttlChunks(i,1) - 1 + bigMod(i) + 1;
            
            devCounter = devCounter + 1;
        end
    else
        disp('No Deviations Found! :D')
    end
    
end

%% Apply TTL QC to spiking data
%now we need to go through the spike data and remove spikes from the
%appropriate periods. 
disp('Applying TTL QC to Spiking Data by Removing Spikes')
for i = 1:numUnits
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    for j = 1:size(devStore,1)
        findSpikes = find(spikeTimes >= dioTimes(devStore(j,4)) & spikeTimes <= dioTimes(devStore(j,5)));
        disp(strcat(num2str(length(findSpikes)),'-spikes removed'))
        spikeTimes(findSpikes) = [];
        s.(desigNames{i}).SpikeTimes = spikeTimes;
    end
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

%% Figure out how DMR file maps onto trodes times
tempVect = linspace(dmrTimes(1),dmrTimes(2),length(stimulus));
for i = 1:length(blockLengths)
    trueDMRtimes = interp1(ttlOnsetTime,dioTimes(ttlChunks(i,1):ttlChunks(i,2)),tempVect);
    dmrStore(:,i) = trueDMRtimes;
end
dmrStep = mean(mode(diff(dmrStore)));
newWin = round(s.Parameters.STRFWin/dmrStep);

%% Now lets start the analysis
%going to go through each thing and perform STRFs

for i = 1:numUnits
    disp(strcat('Working on Unit ',num2str(i)))
    %pulls spike times and times for alignment
    spikeTimes = s.(desigNames{i}).SpikeTimes;
    averageRate = length(spikeTimes) / (s.TimeFilterRange(2) - s.TimeFilterRange(1));
    %make a plot of firing rate over time. 
    sessionFiring = hist(spikeTimes,[s.RotaryData.Velocity(1,1):s.Parameters.SpeedFiringBins:s.RotaryData.Velocity(end,1)]);
    sessionFiring(end) = 0; %this is to compensate for problems with spikes coming after the period
    sessionFiring(1) = 0; %this is to compensate for spikes coming before the tuning period. 
    s.(desigNames{i}).AverageRate = averageRate;
    s.(desigNames{i}).SessionFiring = sessionFiring;
    %lets pull the stuff based on all spikes
    for j = 1:length(blockLengths)
        indivStore = [];
        disp(strcat('Working on DMR Presentation:',num2str(j)))
        %pull out only the spikes during the desired window
        firstSpike = find(spikeTimes > dmrStore(1,j) - s.Parameters.STRFWin(1) ,1,'first');
        lastSpike = find(spikeTimes < dmrStore(end,j) - s.Parameters.STRFWin(2) ,1,'last');
        targetSpikes = spikeTimes(firstSpike:lastSpike);
        for k = 1:length(targetSpikes)
            findTarget = find(dmrStore(:,j) - targetSpikes(k) > 0,1,'first');
            takeChunk = stimulus(:,findTarget+newWin(1):findTarget+newWin(2));
            indivStore(:,:,k) = takeChunk;
        end
        avStore(:,:,i,j) = mean(indivStore,3);
%         bigStore{j} = indivStore;
    end
    bigAvStore(:,:,i) = mean(squeeze(avStore(:,:,i,:)),3);
end



for i = 1:numUnits
    hFig = figure;
    set(hFig, 'Position', [10 80 1240 850])
    %plots average waveform
    subplot(4,4,1)
    hold on
    plot(s.(desigNames{i}).AverageWaveForms,'LineWidth',2)
    title(strcat('AverageFiringRate:',num2str(s.(desigNames{i}).AverageRate)))
    %plots ISI
    subplot(4,4,2)
    hist(s.(desigNames{i}).ISIGraph,1000)
    histMax = max(hist(s.(desigNames{i}).ISIGraph,1000));
    line([s.Parameters.RPVTime s.Parameters.RPVTime],[0 histMax],'LineWidth',1,'Color','red')
    xlim(s.Parameters.ClusterWindow)
    title({strcat('ISI RPV %: ',num2str(s.(desigNames{i}).RPVPercent));...
        strcat(num2str(s.(desigNames{i}).RPVNumber),'/',num2str(s.(desigNames{i}).TotalSpikeNumber))})
    %plot velocity data
    subplot(4,2,2)
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

    subplot(2,2,3)
    imagesc(bigAvStore(:,:,i))
    colormap('parula')
    colorbar
    set(gca,'XTick',[0:60:660]);
    set(gca,'XTickLabel',[-0.1:0.01:0.01]);
    set(gca,'YTick',[1:10:40]);
    set(gca,'YTickLabel',[freqs([1:10:40])]);

    hold off
    spikeGraphName = strcat(fileName,desigNames{i},'SpikeAnalysis');
    savefig(hFig,spikeGraphName);

    %save as PDF with correct name
    set(hFig,'Units','Inches');
    pos = get(hFig,'Position');
    set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hFig,spikeGraphName,'-djpeg','-r0')

end







