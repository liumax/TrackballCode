%This is the Matlab Callback

function toneAllNothingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    stopScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.totalTrials)])))
    scQtUserData.tripSwitch = 1;
    sendScQtControlMessage(['disp(''End of Session'')']);
end


if scQtUserData.trial>=1 && (~isempty(strfind(newLine,'TriggerMatlab')))
    if scQtUserData.trial>1
        scQtUserData.ITI(scQtUserData.trial) = (scQtUserData.lightOn(scQtUserData.trial)...
            -scQtUserData.lightOn(scQtUserData.trial-1))/1000; %This is to log sound off times and ITIs
    end
    
    if ~isfield(scQtUserData,'updateFig') %This code is just to make sure updateFig has a value.
        disp('resetting updateFig');
        scQtUserData.updateFig = -1;
    end
    
    if ~ishandle(scQtUserData.updateFig) %This is to set the basis for all the plots!
        scQtUserData.updateFig = figure('color','w');
        scQtUserData.ax1 = subplot(2,1,1,'parent',scQtUserData.updateFig);
        scQtUserData.ax2 = subplot(2,1,2,'parent',scQtUserData.updateFig);
        hold(scQtUserData.ax1,'on');
        hold(scQtUserData.ax2,'on');
        ylabel(scQtUserData.ax1,'ITI');
        xlabel(scQtUserData.ax1,'Time (s)');
        ylabel(scQtUserData.ax2,'Latency (s)');
        xlabel(scQtUserData.ax2,'Trial Number');
    end

    cla(scQtUserData.ax1); %This clears lickAx. This is because there are constantly adjustments to lickAx
    cla(scQtUserData.ax2);
    %currently can only update the following trial, since figure updates at
    %a time that isnt conducive to this?
    
    %This calculates latency to first lick
    scQtUserData.latTime(scQtUserData.trial) = scQtUserData.stopTime(scQtUserData.trial)...
        - scQtUserData.lightOn(scQtUserData.trial);
    
    title(scQtUserData.ax1,[...
        'Total trials: ', num2str(scQtUserData.trial)]), ...
        
    %plot ITI
    plot(1:scQtUserData.trial,scQtUserData.ITI(1:scQtUserData.trial),'color','k','linewidth',2,...
        'parent',scQtUserData.ax1);
    
    %plot latency
    plot(1:scQtUserData.trial,scQtUserData.latTime(1:scQtUserData.trial),'color','k','linestyle','none',...
        'marker','.','parent',scQtUserData.ax2)
    
    set(scQtUserData.ax1,'ygrid','on');
    set(scQtUserData.ax2,'ygrid','on');

end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage('trigger(3)');
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'Light On'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.lightOn(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1));
end

if (~isempty(strfind(newLine,'Successful Stop'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.stopTime(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1));
end

end