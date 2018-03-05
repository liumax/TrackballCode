%this is testbed for reading out EDR file.

cd F:\170215EDRTests

[data, h] = import_edr('170208_ML170203A_R05_1800_tuningLightStim.EDR');

%this outputs two data: data, which is a nxm array of n inputs and m
%timepoints, with the first column as time, following columns as additional
%information (in this case, 2 is piezo, 3 is TTL)

%h is a structured array with information about settings for the recording
%session. Discard h for now

%parameters
colTime = 1;
colTTL = 3;
colPiezo = 2;

%thresholds
threshTTL = 0.4;

%first, lets pull our TTL data!

edrTTLs = data(:,colTTL);
%find all points above the threshold
ttlFinder = find(edrTTLs>threshTTL);
ttlFinderDiff = diff(ttlFinder);
onsetFinder = [1;find(ttlFinderDiff > 1)+1];

ttlTimes = data(ttlFinder(onsetFinder),colTime); %ttl Times in ms
ttlInds = ttlFinder(onsetFinder);

%lets try making combined average "rasters" based on TTL times
smoothSpan = 5;
sampleRate = 4000; %in Hz
rasterWindow = [-2000 4000]; %lets make this window in samples
rasterHolder = zeros(rasterWindow(2)-rasterWindow(1)+1,length(ttlInds));
for i = 1:length(ttlInds)
    rasterHolder(:,i) = smooth(data(rasterWindow(1)+ttlInds(i):rasterWindow(2)+ttlInds(i),colPiezo),smoothSpan);
end
%downsample!
downSampRaster = downsample(rasterHolder,smoothSpan,2); %the last input is a phase offset. 
rasterVector = [rasterWindow(1):1:rasterWindow(2)]/sampleRate/1000;
rasterVector = downsample(rasterVector,smoothSpan,2); 

figure
hold on
for i = 1:length(ttlInds)
    plot(rasterHolder(:,i)+(i*0.005),'Color',rand(3,1))
end
plot([2000 2000],[0 800*0.01])
plot([2400 2400],[0 800*0.01])
%this shows a ton of responses. A fair number of the tones appear to
%produce responses

absHolder = abs(rasterHolder);
figure
plot(mean(absHolder,2));
%absolute value shows some kind of response post tone!!

%i think moving forward, the best plan of action is to smooth over 5
%samples (basically averaging over 1ms) and then downsample by factor of 4
%(using matlab function DOWNSAMPLE). This will save on memory use, while
%still providing a useful way to align things. 

%I think I should do this with the data initially? actually that doesnt
%really make sense...

% OK lets do this: we'll do everything as normal. Then we downsample the
% rasters and also the overall traces, saving the middle time? or something
% like that...

rasterHolder = smooth(rasterHolder,5);


