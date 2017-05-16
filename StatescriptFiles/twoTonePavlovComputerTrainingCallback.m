%This is the Matlab Callback

function twoTonePavlovComputerTrainingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    twoTonePavlovComputerScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

%this kills things after the last trial of the session.
if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.totalTrials)])))
    scQtUserData.tripSwitch = 1;
end

%this is for the first trial, which is triggered via the script. This sends
%information to the mbed for execution. THIS DOES NOT PLAY THE TONE
if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.Master(scQtUserData.trial,1))]); 
    sendScQtControlMessage('trigger(1)');
end
%This is for all other trials!
if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.Master(scQtUserData.trial,1))]); 
    sendScQtControlMessage('trigger(1)');
end

if ~isempty(strfind(newLine,'TriggerSound'))
    soundID = scQtUserData.Master(scQtUserData.trial,2);
    if sound ID == 1
        sound(scQtUserData.ToneBig,192000)
    elseif sound ID == 2
        sound(scQtUserData.ToneSmall,192000)
    end
end

end