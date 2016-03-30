%This is the Matlab Callback

function lickTrainingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    lickTrainingScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.blocks*scQtUserData.blockSize)])))
    scQtUserData.tripSwitch = 1;
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.master(scQtUserData.trial,2))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.master(scQtUserData.trial,3))]);
    sendScQtControlMessage(['rewLength = ',num2str(scQtUserData.master(scQtUserData.trial,1))]);
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.master(scQtUserData.trial,2))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.master(scQtUserData.trial,3))]);
    sendScQtControlMessage(['rewLength = ',num2str(scQtUserData.master(scQtUserData.trial,1))]);
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'SoundOn'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.master(scQtUserData.trial,4) = str2num(newLine(1:spaceFinder(1)-1));
end

end