%This is the Matlab Callback

function twoTonePavlovComputerTrainingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    pavThreeToneWAVScript;
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
    try
        scQtUserData.licks(scQtUserData.lickCounter,1) = str2num(newLine(1:spaceFinder(1)-1))/1000;
        scQtUserData.licks(scQtUserData.lickCounter,2) = scQtUserData.trial;
        scQtUserData.licks(scQtUserData.lickCounter,4) = scQtUserData.LickDesig;
        scQtUserData.lickCounter = scQtUserData.lickCounter + 1;
    catch
    end
end

if(~isempty(strfind(newLine,'Tone Delivered')))
    spaceFinder = find(newLine == ' ');
    scQtUserData.cueTime(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1))/1000;
end

%this is for the first trial, which is triggered via the script. This sends
%information to the mbed for execution. THIS DOES NOT PLAY THE TONE
if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(round(scQtUserData.Master(scQtUserData.trial,1)))]);
    sendScQtControlMessage(['outDur = ',num2str(round(scQtUserData.Master(scQtUserData.trial,3)))]);  
    sendScQtControlMessage(['toneOutDel =',num2str(scQtUserData.outDelayMatrix(scQtUserData.trial))]);
    sendScQtControlMessage(['outPort = ',num2str(round(scQtUserData.Master(scQtUserData.trial,4)))]);
    sendScQtControlMessage(['soundPort = ',num2str(round(scQtUserData.Master(scQtUserData.trial,5)))]);  
    scQtUserData.LickDesig = scQtUserData.Master(scQtUserData.trial,2);
    soundID = scQtUserData.Master(scQtUserData.trial,2)
    if soundID == 1
        sendScQtControlMessage(['disp(''PlaySmall'')']);
    elseif soundID == 2
        sendScQtControlMessage(['disp(''PlayBig'')']);
    elseif soundID == 3
        sendScQtControlMessage(['disp(''Play Pun'')']);
    elseif soundID == 4
        sendScQtControlMessage(['disp(''PlayFreeRew'')']);
    elseif soundID == 5
        sendScQtControlMessage(['disp(''PlayBigCatch'')']);
    elseif soundID == 6
        sendScQtControlMessage(['disp(''PlayFreePun'')']);
    elseif soundID == 7
        sendScQtControlMessage(['disp(''PlayPunCatch'')']);
    end
    sendScQtControlMessage('trigger(1)');
end
%This is for all other trials!

if ~isempty(strfind(newLine,'PlotTime')) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(round(scQtUserData.Master(scQtUserData.trial,1)))]);
    sendScQtControlMessage(['outDur = ',num2str(round(scQtUserData.Master(scQtUserData.trial,3)))]);  
    sendScQtControlMessage(['toneOutDel =',num2str(scQtUserData.outDelayMatrix(scQtUserData.trial))]);
    sendScQtControlMessage(['outPort = ',num2str(round(scQtUserData.Master(scQtUserData.trial,4)))]);
    sendScQtControlMessage(['soundPort = ',num2str(round(scQtUserData.Master(scQtUserData.trial,5)))]);  
    scQtUserData.LickDesig = scQtUserData.Master(scQtUserData.trial,2);
    soundID = scQtUserData.Master(scQtUserData.trial,2)
    if soundID == 1
        sendScQtControlMessage(['disp(''PlaySmall'')']);
    elseif soundID == 2
        sendScQtControlMessage(['disp(''PlayBig'')']);
    elseif soundID == 3
        sendScQtControlMessage(['disp(''Play Pun'')']);
    elseif soundID == 4
        sendScQtControlMessage(['disp(''PlayFreeRew'')']);
    elseif soundID == 5
        sendScQtControlMessage(['disp(''PlayBigCatch'')']);
    elseif soundID == 6
        sendScQtControlMessage(['disp(''PlayFreePun'')']);
    elseif soundID == 7
        sendScQtControlMessage(['disp(''PlayPunCatch'')']);
    end
    sendScQtControlMessage('trigger(1)');
    
    try
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


        %cleans up raster with the cue time
        scQtUserData.licks(scQtUserData.licks(:,2) == scQtUserData.trial,1) = scQtUserData.licks(scQtUserData.licks(:,2) == scQtUserData.trial,1) - scQtUserData.cueTime(scQtUserData.trial);

        %plot things!
        subplot(3,1,1)
        lickTrunc = scQtUserData.licks(1:scQtUserData.lickCounter,:);
%         findLow = find(lickTrunc(:,4) == 1);
%         findHi = find(lickTrunc(:,4) == 2);
        
        plot(lickTrunc(:,1),lickTrunc(:,2),'b.',...
            'parent',scQtUserData.ax1);
        axis(scQtUserData.ax1,[-2 6 0.5 200])
        title(num2str(scQtUserData.trial))

%         plot(lickTrunc(findLow,1),lickTrunc(findLow,2),'b.',...
%             'parent',scQtUserData.ax1);
%         hold on
%         plot(lickTrunc(findHi,1),lickTrunc(findHi,2),'r.',...
%             'parent',scQtUserData.ax1);
%         axis(scQtUserData.ax1,[scQtUserData.lickAxes(1) scQtUserData.lickAxes(end) 0.5 scQtUserData.trial])

        subplot(3,1,2)

        subplot(3,1,3)
        %process histogram data
    %     histData = scQtUserData.licks(1:scQtUserData.lickCounter,1);
        lowLicks = lickTrunc(lickTrunc(:,4) == 1,1);
        lowLicks(lowLicks > scQtUserData.lickAxes(end)) = [];
        lowLicks(lowLicks < scQtUserData.lickAxes(1)) = [];
        lowHist = hist(lowLicks,scQtUserData.lickAxes);

        hiLicks = lickTrunc(lickTrunc(:,4) == 2,1);
        hiLicks(hiLicks > scQtUserData.lickAxes(end)) = [];
        hiLicks(hiLicks < scQtUserData.lickAxes(1)) = [];
        hiHist = hist(hiLicks,scQtUserData.lickAxes);
        
        freeLicks = lickTrunc(lickTrunc(:,4) == 3,1);
        freeLicks(freeLicks > scQtUserData.lickAxes(end)) = [];
        freeLicks(freeLicks < scQtUserData.lickAxes(1)) = [];
        freeHist = hist(freeLicks,scQtUserData.lickAxes);
        
        catchLicks = lickTrunc(lickTrunc(:,4) == 4,1);
        catchLicks(catchLicks > scQtUserData.lickAxes(end)) = [];
        catchLicks(catchLicks < scQtUserData.lickAxes(1)) = [];
        catchHist = hist(catchLicks,scQtUserData.lickAxes);

        hold on
        plot(scQtUserData.lickAxes,lowHist,'r')
        plot(scQtUserData.lickAxes,hiHist,'b')
        plot(scQtUserData.lickAxes,freeHist,'g.-')
        plot(scQtUserData.lickAxes,catchHist,'k.-')
    catch
        subplot(3,1,1)
        title(num2str(scQtUserData.trial))
        disp(strcat('lengthLicks',num2str(length(scQtUserData.licks))));
        disp(strcat('lickCounter',num2str(scQtUserData.lickCounter)));
    end
    
end

end