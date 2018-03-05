%This is the testbed for code to examine basic one tone licking behavior

fname = '170607ML170504H';
rasterWindow = [-10 10];

%First thing I need to do is extract the file! 
[trialStates, portStates, trialParams] = maxTrialVariablesLickingTask(fname); 

%next is to extract rotary data
[x] = functionMBEDrotary(portStates.inStates(:,4),portStates.inStates(:,5),portStates.tStamps/1000,0.1);

%now lets pull tone delivery times
cueTimes = trialStates.cueTimeRew1/1000;

%pull lick times
lickTimes = trialParams.licking(:,1)/1000;

%make licking raster
[lickRaster] = functionBasicRaster(lickTimes,cueTimes',rasterWindow);
figure
plot(lickRaster(:,1),lickRaster(:,2),'b.')


%now make raster of velocity
velAxis = [rasterWindow(1):0.1:rasterWindow(2)];
velRaster = zeros(length(velAxis),length(cueTimes));
for i = 1:length(cueTimes)
    %first find the nearest time point
    findTime = find(x.Velocity(:,1) - cueTimes(i) > 0,1,'first');
    velRaster(:,i) = x.Velocity(findTime:findTime+length(velAxis)-1,2);
end






