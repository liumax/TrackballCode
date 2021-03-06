%This is the Matlab Callback

function oneTonePavlovComputerCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    oneTonePavlovComputerScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

%this kills things after the last trial of the session.
if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.totalTrials)])))
    scQtUserData.tripSwitch = 1;
    sendScQtControlMessage(['disp(''EndSession'')']);
end

if(~isempty(strfind(newLine,'Lick Detected')))
    spaceFinder = find(newLine == ' ');
    scQtUserData.licks(scQtUserData.lickCounter,1) = str2num(newLine(1:spaceFinder(1)-1))/1000;
    scQtUserData.licks(scQtUserData.lickCounter,2) = scQtUserData.trial;
%     scQtUserData.licks(scQtUserData.lickCounter,3) = scQtUserData.trPhase;
    scQtUserData.lickCounter = scQtUserData.lickCounter + 1;
end

if(~isempty(strfind(newLine,'Tone Delivered')))
    spaceFinder = find(newLine == ' ');
    scQtUserData.cueTime(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1))/1000;
end

if(~isempty(strfind(newLine,'trPhase')))
    spaceFinder = find(newLine == ' ');
    scQtUserData.trPhase = str2num(newLine(spaceFinder(2)+1:end));
end

%this is for the first trial, which is triggered via the script. This sends
%information to the mbed for execution. THIS DOES NOT PLAY THE TONE
if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(round(scQtUserData.Master(scQtUserData.trial,1)))]);
    sendScQtControlMessage(['rewLength = ',num2str(round(scQtUserData.Master(scQtUserData.trial,3)))]);  
    sendScQtControlMessage('trigger(1)');
end
%This is for all other trials!
if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(round(scQtUserData.Master(scQtUserData.trial,1)))]); 
    sendScQtControlMessage(['rewLength = ',num2str(round(scQtUserData.Master(scQtUserData.trial,3)))]);  
    sendScQtControlMessage('trigger(1)');
end

if ~isempty(strfind(newLine,'TriggerSound'))
    disp('PlayBig')
    sound(scQtUserData.ToneBig,192000)
end

if ~isempty(strfind(newLine,'PlotTime'))
    if ~isfield(scQtUserData,'updateFig') %This code is just to make sure updateFig has a value.
        disp('resetting updateFig');
        scQtUserData.updateFig = -1;
    end
    
    
    if ~ishandle(scQtUserData.updateFig) %This is to set the basis for all the plots!
        scQtUserData.updateFig = figure('color','w');
        scQtUserData.ax1 = subplot(3,1,1,'parent',scQtUserData.updateFig);
        scQtUserData.ax2 = subplot(3,1,2,'parent',scQtUserData.updateFig);
        scQtUserData.ax3 = subplot(3,1,3,'parent',scQtUserData.updateFig);
        hold(scQtUserData.ax1,'on');
        hold(scQtUserData.ax2,'on');
        hold(scQtUserData.ax3,'on');
        ylabel(scQtUserData.ax1,'Trial #');
        xlabel(scQtUserData.ax1,'Time (s)');
        ylabel(scQtUserData.ax2,'Reaction Time (s)');
        xlabel(scQtUserData.ax2,'Time (s)');
        ylabel(scQtUserData.ax3,'Licks');
        xlabel(scQtUserData.ax3,'Time (s)');
    end
    
    cla(scQtUserData.ax1);
    cla(scQtUserData.ax2); %clears ax2
    cla(scQtUserData.ax3);
    
    %finds first lick and stores
    firstLick = scQtUserData.licks(:,1) - scQtUserData.cueTime(scQtUserData.trial);
    firstLick(firstLick <= 0) = [];
    if numel(firstLick)~= 0
        scQtUserData.firstLick(scQtUserData.trial) = firstLick(1);
    else
        scQtUserData.firstLick(scQtUserData.trial) = 0;
    end
    
    %cleans up raster with the cue time
    scQtUserData.licks(scQtUserData.licks(:,2) == scQtUserData.trial,1) = scQtUserData.licks(scQtUserData.licks(:,2) == scQtUserData.trial,1) - scQtUserData.cueTime(scQtUserData.trial);
    
    %plot things!
    subplot(3,1,1)
    plot(scQtUserData.licks(1:scQtUserData.lickCounter,1),scQtUserData.licks(1:scQtUserData.lickCounter,2),'b.',...
        'parent',scQtUserData.ax1);
    axis(scQtUserData.ax1,[scQtUserData.lickAxes(1) scQtUserData.lickAxes(end) 0.5 scQtUserData.trial])
    
    subplot(3,1,2)
    plot(scQtUserData.firstLick,'b.','parent',scQtUserData.ax2)
    
    subplot(3,1,3)
    %process histogram data
    histData = scQtUserData.licks(1:scQtUserData.lickCounter,1);
    histData(histData > scQtUserData.lickAxes(end)) = [];
    histData(histData < scQtUserData.lickAxes(1)) = [];
    hist(histData,scQtUserData.lickAxes)
    
end

end