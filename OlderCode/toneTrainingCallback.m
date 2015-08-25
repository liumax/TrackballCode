%This is the Matlab Callback

function toneTrainingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    toneTrainingScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.totalTrials)])))
    scQtUserData.tripSwitch = 1;
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.ITI(scQtUserData.trial))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.soundRewDel(scQtUserData.trial))]);
    sendScQtControlMessage(['sound = ',num2str(scQtUserData.soundDeliver(scQtUserData.trial))]);
    if scQtUserData.toneRule == 1
        sendScQtControlMessage('trigger(1)');
    elseif scQtUserData.toneRule == 0
        sendScQtControlMessage('trigger(3)');
    end
    
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.ITI(scQtUserData.trial))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.soundRewDel(scQtUserData.trial))]);
    sendScQtControlMessage(['sound = ',num2str(scQtUserData.soundDeliver(scQtUserData.trial))]);
    if scQtUserData.toneRule == 1
        sendScQtControlMessage('trigger(1)');
    elseif scQtUserData.toneRule == 0
        sendScQtControlMessage('trigger(3)');
    end
end

end