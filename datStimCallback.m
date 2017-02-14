%This is the Matlab Callback

function tuningCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    minITI = 10; %minimum ITI in seconds
    maxITI = 15; %maximum ITI in seconds
    k = 2.5;
    p = (1-exp(-k))*rand(1000,1);
    tau = (maxITI-minITI)/k;
    x = round(minITI + (-log(1-p))*tau); 
    scQtUserData.ITIs = x*1000;
    scQtInitiated = 1;
    newLine = 'start next trial';
end


if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['ITI = ',num2str(scQtUserData.ITIs(scQtUserData.trial))]); 
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    toneWave = sin(2*pi*(toneFreq/scQtUserData.fs)*(1:scQtUserData.L))';
    finalWave = toneWave.*scQtUserData.rampProfile;
    finalWave = finalWave*toneAmpl;
    paddedWave = zeros(scQtUserData.paddingL,1);
    paddedWave(1:size(finalWave,1)) = finalWave;
    scQtUserData.soundVector = [paddedWave,scQtUserData.ttlSig];
    
    sendScQtControlMessage(['trigger(3)']);
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.Master(scQtUserData.trial,4))]); 
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    toneFreq = scQtUserData.Master(scQtUserData.trial,1);
    toneAmpl = scQtUserData.Master(scQtUserData.trial,3);
    toneWave = sin(2*pi*(toneFreq/scQtUserData.fs)*(1:scQtUserData.L))';
    finalWave = toneWave.*scQtUserData.rampProfile;
    finalWave = finalWave*toneAmpl;
    paddedWave = zeros(scQtUserData.paddingL,1);
    paddedWave(1:size(finalWave,1)) = finalWave;
    scQtUserData.soundVector = [paddedWave,scQtUserData.ttlSig];
    
    sendScQtControlMessage(['trigger(3)']);
end

if (~isempty(strfind(newLine,'Play Sound')))
    sound(scQtUserData.soundVector,scQtUserData.fs);
    sendScQtControlMessage(['disp(''Trial: ',num2str(scQtUserData.trial),'/',num2str(length(scQtUserData.Master)),' Frequency:',num2str(scQtUserData.Master(scQtUserData.trial,1)),' DB:',num2str(scQtUserData.Master(scQtUserData.trial,2)),''')']);
end

end