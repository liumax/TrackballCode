

%This is the Matlab Callback

function yangDanOperantTrainingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    audLickOperantTrainingScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.totalTrials)])))
    scQtUserData.tripSwitch = 1;
    sendScQtControlMessage(['disp(''End of Session'')']);
end

if scQtUserData.trial>=1 && (~isempty(strfind(newLine,'TriggerMatlab')))
    
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
        hold(scQtUserData.ax2,'on');
        ylabel(scQtUserData.ax1,'Trial #');
        xlabel(scQtUserData.ax1,'Time (s)');
        ylabel(scQtUserData.ax2,'Cumulative Licks');
        xlabel(scQtUserData.ax2,'Time (s)');
        ylabel(scQtUserData.ax3,'Reaction Time (s)');
        xlabel(scQtUserData.ax3,'Trial Number');
    end
    
    cla(scQtUserData.ax1);
    cla(scQtUserData.ax2); %clears ax2
    
    %this preps licking data by effectively binning it
       scQtUserData.LickTime = round((scQtUserData.LickTime - scQtUserData.cueTime(scQtUserData.trial))/scQtUserData.binSize)...
        +((scQtUserData.histLim(1)*-1)*1000/scQtUserData.binSize); %note the fudge factor here
    scQtUserData.LickTime(scQtUserData.LickTime > length(scQtUserData.graphAxes)) = [];
    scQtUserData.LickTime(scQtUserData.LickTime < 1) = [];
    
    %this sorts data into respective rasters and histograms
    scQtUserData.rewHist([scQtUserData.LickTime],1)=scQtUserData.rewHist([scQtUserData.LickTime],1)+1;
    scQtUserData.rewRast(scQtUserData.rewHolder:scQtUserData.rewHolder+length(scQtUserData.LickTime)-1,1) = scQtUserData.trial;
    scQtUserData.rewRast(scQtUserData.rewHolder:scQtUserData.rewHolder+length(scQtUserData.LickTime)-1,2) = (scQtUserData.LickTime-20)/1000*scQtUserData.binSize;
    scQtUserData.rewHolder = scQtUserData.rewHolder + length(scQtUserData.LickTime);
    
    title(scQtUserData.ax1,[...
        'Total trials: ', num2str(scQtUserData.trial)]), ...

    %plot raster
    plot(scQtUserData.rewRast(:,2),scQtUserData.rewRast(:,1),'g.',...
        'parent',scQtUserData.ax1);
    axis(scQtUserData.ax1,[scQtUserData.histLim(1) scQtUserData.histLim(2) 0.5 scQtUserData.trial])

    %plot histogram
    plot(scQtUserData.graphAxes,scQtUserData.rewHist,'color','g','linewidth',1,...
        'parent',scQtUserData.ax2);
    
    %plot reaction timing
    plot(1:1:scQtUserData.totalTrials,scQtUserData.rewTime,'k.',...
        'parent',scQtUserData.ax3);
    
    set(scQtUserData.ax1,'ygrid','on');
    set(scQtUserData.ax2,'ygrid','on');
    set(scQtUserData.ax3,'ygrid','on');
    
    scQtUserData.LickTime = zeros (1000,1);
    scQtUserData.lickHolder = 1;
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage('trigger(2)');
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'Cue On')));
    spaceFinder = find(newLine == ' ');
    scQtUserData.cueTime(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1));
end

if (~isempty(strfind(newLine,'Lick Detected')));
    spaceFinder = find(newLine == ' ');
    scQtUserData.LickTime(scQtUserData.lickHolder) = str2num(newLine(1:spaceFinder(1)-1));
    scQtUserData.lickHolder = scQtUserData.lickHolder + 1;
end

if (~isempty(strfind(newLine,'Reward Delivered')));
    spaceFinder = find(newLine == ' ');
    scQtUserData.rewTime(scQtUserData.trial) = (str2num(newLine(1:spaceFinder(1)-1))...
        -scQtUserData.cueTime(scQtUserData.trial))/1000;
end

end