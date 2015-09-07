%This is the Matlab Callback

function yangDan1Callback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    yangDan1Script;
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
        scQtUserData.ax1 = subplot(2,1,1,'parent',scQtUserData.updateFig);
        scQtUserData.ax2 = subplot(2,1,2,'parent',scQtUserData.updateFig);
        hold(scQtUserData.ax1,'on');
        hold(scQtUserData.ax2,'on');
        ylabel(scQtUserData.ax1,'Trial #');
        xlabel(scQtUserData.ax1,'Time (s)');
        ylabel(scQtUserData.ax2,'Licks/Sec');
        xlabel(scQtUserData.ax2,'Time (s)');
    end
    
    cla(scQtUserData.ax1);
    cla(scQtUserData.ax2); %clears ax2
    
    %this preps licking data by effectively binning it
    scQtUserData.LickTime = round((scQtUserData.LickTime - scQtUserData.cueTime(scQtUserData.trial))/scQtUserData.binSize)+20; %note the 20 fudge factor here
    scQtUserData.LickTime(scQtUserData.LickTime > 70) = [];
    scQtUserData.LickTime(scQtUserData.LickTime < 1) = [];
    
    %this sorts data into respective rasters and histograms
    if scQtUserData.sessionType(scQtUserData.trial) == 1
        scQtUserData.neutHist([scQtUserData.LickTime],1)=scQtUserData.neutHist([scQtUserData.LickTime],1)+1;
        scQtUserData.neutRast(scQtUserData.neutHolder:scQtUserData.neutHolder+length(scQtUserData.LickTime)-1,1) = scQtUserData.trial;
        scQtUserData.neutRast(scQtUserData.neutHolder:scQtUserData.neutHolder+length(scQtUserData.LickTime)-1,2) = (scQtUserData.LickTime-20)/1000*scQtUserData.binSize;
        scQtUserData.neutHolder = scQtUserData.neutHolder + length(scQtUserData.LickTime);
    elseif scQtUserData.sessionType(scQtUserData.trial) == 2
        scQtUserData.rewHist([scQtUserData.LickTime],1)=scQtUserData.rewHist([scQtUserData.LickTime],1)+1;
        scQtUserData.rewRast(scQtUserData.rewHolder:scQtUserData.rewHolder+length(scQtUserData.LickTime)-1,1) = scQtUserData.trial;
        scQtUserData.rewRast(scQtUserData.rewHolder:scQtUserData.rewHolder+length(scQtUserData.LickTime)-1,2) = (scQtUserData.LickTime-20)/1000*scQtUserData.binSize;
        scQtUserData.rewHolder = scQtUserData.rewHolder + length(scQtUserData.LickTime);
    elseif scQtUserData.sessionType(scQtUserData.trial) == 3
        scQtUserData.punHist([scQtUserData.LickTime],1)=scQtUserData.punHist([scQtUserData.LickTime],1)+1;
        scQtUserData.punRast(scQtUserData.punHolder:scQtUserData.punHolder+length(scQtUserData.LickTime)-1,1) = scQtUserData.trial;
        scQtUserData.punRast(scQtUserData.punHolder:scQtUserData.punHolder+length(scQtUserData.LickTime)-1,2) = (scQtUserData.LickTime-20)/1000*scQtUserData.binSize;
        scQtUserData.punHolder = scQtUserData.punHolder + length(scQtUserData.LickTime);
    end
    
    title(scQtUserData.ax1,[...
        'Total trials: ', num2str(scQtUserData.trial)]), ...

    %plot raster
    plot(scQtUserData.neutRast(:,2),scQtUserData.neutRast(:,1),'b.',...
        'parent',scQtUserData.ax1);
    plot(scQtUserData.rewRast(:,2),scQtUserData.rewRast(:,1),'g.',...
        'parent',scQtUserData.ax1);
    plot(scQtUserData.punRast(:,2),scQtUserData.punRast(:,1),'r.',...
        'parent',scQtUserData.ax1);
    plot(scQtUserData.rewTime(:,2),scQtUserData.rewTime(:,1),'go',...
        'parent',scQtUserData.ax1);
    plot(scQtUserData.punTime(:,2),scQtUserData.punTime(:,1),'ro',...
        'parent',scQtUserData.ax1);
    axis(scQtUserData.ax1,[-2 5 0.5 scQtUserData.trial])
    
    %plot histogram
     plot(scQtUserData.graphAxes,scQtUserData.neutHist,'color','b','linewidth',1,...
        'parent',scQtUserData.ax2);
     plot(scQtUserData.graphAxes,scQtUserData.rewHist,'color','g','linewidth',1,...
        'parent',scQtUserData.ax2);
    plot(scQtUserData.graphAxes,scQtUserData.punHist,'color','r','linewidth',1,...
        'parent',scQtUserData.ax2);
    
    set(scQtUserData.ax1,'ygrid','on');
    set(scQtUserData.ax2,'ygrid','on');
    
    scQtUserData.LickTime = zeros (1000,1);
    scQtUserData.lickHolder = 1;
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage(['trigger(',num2str(scQtUserData.sessionType(scQtUserData.trial)),')']);
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.itiTime(scQtUserData.trial))]); 
    sendScQtControlMessage(['trigger(',num2str(scQtUserData.sessionType(scQtUserData.trial)),')']);
end

if (~isempty(strfind(newLine,'Cue Light On')));
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
    scQtUserData.rewTime(scQtUserData.trial,1) = scQtUserData.trial;
    scQtUserData.rewTime(scQtUserData.trial,2) = round((str2num(newLine(1:spaceFinder(1)-1))-...
        scQtUserData.cueTime(scQtUserData.trial))/scQtUserData.binSize)/1000*scQtUserData.binSize;
end

if (~isempty(strfind(newLine,'Punishment Delivered')));
    spaceFinder = find(newLine == ' ');
    scQtUserData.punTime(scQtUserData.trial,1) = scQtUserData.trial;
    scQtUserData.punTime(scQtUserData.trial,2) = round((str2num(newLine(1:spaceFinder(1)-1))-...
        scQtUserData.cueTime(scQtUserData.trial))/scQtUserData.binSize)/1000*scQtUserData.binSize;
end

end