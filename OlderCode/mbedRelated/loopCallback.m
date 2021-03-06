function loopCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    loopScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

pause(10)
sendScQtControlMessage('trigger(1)');
sendScQtControlMessage('trigger(2)');