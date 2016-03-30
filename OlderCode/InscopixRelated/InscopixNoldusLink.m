%%%%%%%%%%%%%%%%%%%%%%%%TrackBall Code%%%%%%%%%%%%%%%%%%%%%%%%

recordtime=22; %total amount of time recording in minutes
recordblock=5; %amount of time per block in seconds
numblocks=recordtime*60/recordblock;
master=zeros(numblocks,1);
master(:,1)=recordblock;
master(:,1)=cumsum(master(:,1));

% totaltime=35; %total time to get samples for in minutes
% totaltime=totaltime*60;
totaltime=recordtime*60;

[fname pname] = uiputfile('test1.mat');

% ao= analogoutput('nidaq','Dev1');
% ch= addchannel(ao,0,'Inscopix');

ai = analoginput('nidaq','Dev3');
inscopix=addchannel(ai,0,'InscopixFeed');
noldus=addchannel(ai,1,'Noldus');

sampleRate = get(ai, 'SampleRate')
ai.SamplesPerTrigger = sampleRate*(totaltime);
get(ai, 'SamplesPerTrigger')
    
InscopixData=zeros(sampleRate*(totaltime),3);
InscopixTime=zeros(sampleRate*(totaltime),1);

t0=clock;

start(ai);
currPtAi=1;

% putsample(ao,5);
    figure
for q=1:numblocks

%     rnow=etime(clock,t0);
%     del=master(numblocks,1)-rnow;
%     pause(del);
    pause(recordblock);
    aiSamples = ai.SamplesAvailable
    [InscopixData((currPtAi):(currPtAi-1+aiSamples),1:2), InscopixTime((currPtAi):(currPtAi-1+aiSamples))] = getdata(ai,aiSamples);
    clf;                                    %Clears figure
    plot(InscopixTime((currPtAi):(currPtAi-1+aiSamples)),InscopixData((currPtAi):(currPtAi-1+aiSamples),1:2)),grid on         %Updates figure with new data
    title('All Data')
    currPtAi=currPtAi+aiSamples;            %Updates currPt counter
end

% putsample(ao,0);
aiSamples = ai.SamplesAvailable
[InscopixData((currPtAi):(currPtAi-1+aiSamples),1:2), InscopixTime((currPtAi):(currPtAi-1+aiSamples))] = getdata(ai,aiSamples);

stop(ai);
disp 'End of Recording Period'

save(fullfile(pname,fname),'InscopixTime','InscopixData','recordtime','sampleRate')

