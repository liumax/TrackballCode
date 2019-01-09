

%% open the s file



%% navigate to LFP folder

lfpWindow = [-0.2 0.5];

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

%for plotting different electrode LFPs over all tones
for i = 1:length(lfpFiles)
    [lfpTrace lfpTimes lfpSampleRate] = functionLFPExtractor((lfpFiles{i}));
    [lfpAverage,lfpSTD] = functionLFPpsth(lfpTimes,lfpTrace,lfpWindow,lfpSampleRate,s.TTLs.Tuning1);
    figure
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],lfpAverage)
    hold on
    [lfpAverage,lfpSTD] = functionLFPpsth(lfpTimes,lfpTrace,lfpWindow,lfpSampleRate,s.TTLs.Tuning2);
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],lfpAverage,'r')
%     hold on
%     plot([-1:1/lfpSampleRate:3-1/lfpSampleRate],lfpAverage-lfpSTD/sqrt(200))
%     plot([-1:1/lfpSampleRate:3-1/lfpSampleRate],lfpAverage+lfpSTD/sqrt(200))
end

%get the indices of the first tuning by frequency

tuningFreqs = s.SoundData.Tuning1.UniqueFreqs;
tuning1 = zeros(length(s.SoundData.Tuning1.Frequencies)/length(tuningFreqs),length(tuningFreqs));
tuning2 = zeros(length(s.SoundData.Tuning2.Frequencies)/length(tuningFreqs),length(tuningFreqs));

for i = 1:length(tuningFreqs)
    tuning1(:,i) = find(s.SoundData.Tuning1.Frequencies == tuningFreqs(i));
    tuning2(:,i) = find(s.SoundData.Tuning2.Frequencies == tuningFreqs(i));
end

tuning1LFPHolder = cell(length(tuningFreqs),numLFPs);
tuning2LFPHolder = cell(length(tuningFreqs),numLFPs);
for i = 1:numLFPs
    [lfpTrace lfpTimes lfpSampleRate] = functionLFPExtractor((lfpFiles{i}));
    for j = 1:length(tuningFreqs)
        [lfpAverage,lfpSTD] = functionLFPpsth(lfpTimes,lfpTrace,lfpWindow,lfpSampleRate,s.TTLs.Tuning1(tuning1(:,j)));
        tuning1LFPHolder{j,i} = lfpAverage;
        [lfpAverage,lfpSTD] = functionLFPpsth(lfpTimes,lfpTrace,lfpWindow,lfpSampleRate,s.TTLs.Tuning2(tuning2(:,j)));
        tuning2LFPHolder{j,i} = lfpAverage;
    end
end

for i = 1:numLFPs
    figure
    subplot(4,1,1)
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning1LFPHolder{1,i})
    hold on
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning2LFPHolder{1,i},'r')
    title('4kHz')
    
    subplot(4,1,2)
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning1LFPHolder{2,i})
    hold on
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning2LFPHolder{2,i},'r')
    title('8kHz')
    
    subplot(4,1,3)
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning1LFPHolder{3,i})
    hold on
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning2LFPHolder{3,i},'r')
    title('16kHz')
    
    subplot(4,1,4   )
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning1LFPHolder{4,i})
    hold on
    plot([lfpWindow(1):1/lfpSampleRate:lfpWindow(2)-1/lfpSampleRate],tuning2LFPHolder{4,i},'r')
    title('32kHz')
end

