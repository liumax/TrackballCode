%This is the Matlab Callback

function recordToneCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    recordToneScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,'Recording Trials Begin Now')))
    scQtUserData.triggerSwitch = 1
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.triggerSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.triggerSwitch == 1;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]);
    if scQtUserData.triggerCounter < scQtUserData.freqStim
        scQtUserData.triggerCounter = scQtUserData.triggerCounter + 1;
        sendScQtControlMessage('trigger(1)');
    elseif scQtUserData.triggerCounter >= scQtUserData.freqStim
        scQtUserData.triggerCounter = 1;
        sendScQtControlMessage('trigger(2)');
    end
end

end