%This is the Matlab Callback

function tuningCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    tuningScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(length(scQtUserData.Master))])))
    scQtUserData.tripSwitch = 1;
    sendScQtControlMessage(['disp(''End of Session'')']);
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.Master(scQtUserData.trial,4))]); 
    sendScQtControlMessage(['trigger(3)']);
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.Master(scQtUserData.trial,4))]); 
    sendScQtControlMessage(['trigger(3)']);
end

if (~isempty(strfind(newLine,'Play Sound')))
    toneFreq = scQtUserData.Master(scQtUserData.trial,1);
    toneAmpl = scQtUserData.Master(scQtUserData.trial,3);
    toneWave = sin(2*pi*(toneFreq/scQtUserData.fs)*(1:scQtUserData.L))';
    finalWave = toneWave.*scQtUserData.rampProfile;
    finalWave = finalWave*toneAmpl;
    paddedWave = zeros(scQtUserData.paddingL,1);
    paddedWave(1:size(finalWave,1)) = finalWave;
    soundVector = [paddedWave,scQtUserData.ttlSig];
    sound(soundVector,scQtUserData.fs);
    sendScQtControlMessage(['disp(''Trial: ',num2str(scQtUserData.trial),'/',num2str(length(scQtUserData.Master)),' Frequency:',num2str(scQtUserData.Master(scQtUserData.trial,1)),' DB:',num2str(scQtUserData.Master(scQtUserData.trial,2)),''')']);
end

end