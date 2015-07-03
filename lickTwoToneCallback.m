%This is the Matlab Callback

function lickTrainingCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    lickTwoToneScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.totalTrials)])))
    scQtUserData.tripSwitch = 1;
end

if ~isempty(strfind(newLine,'Lick Detected')) && scQtUserData.tripSwitch == 0;
    spaceFinder= find(newLine == ' ');
    scQtUserData.licks(scQtUserData.lickCounter,1) = str2num(newLine(1:spaceFinder(1)-1));
    scQtUserData.lickCounter= scQtUserData.lickCounter+1;
end

if scQtUserData.trial>=1 && (~isempty(strfind(newLine,'TriggerMatlab')))
    if scQtUserData.trial>1
        scQtUserData.ITI(scQtUserData.trial) = (scQtUserData.soundOn(scQtUserData.trial)...
            -scQtUserData.soundOn(scQtUserData.trial-1))/1000; %This is to log sound off times and ITIs
    end
    
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
        ylabel(scQtUserData.ax1,'Licks');
        xlabel(scQtUserData.ax1,'Time from Cue Onset (s)');
        ylabel(scQtUserData.ax2,'Latency (s)');
        xlabel(scQtUserData.ax2,'Trial Number');
        ylabel(scQtUserData.ax3,'Anticipatory Licks');
        xlabel(scQtUserData.ax3,'Trial Number');
    end

    cla(scQtUserData.ax1); %This clears lickAx. This is because there are constantly adjustments to lickAx
    cla(scQtUserData.ax2);
    cla(scQtUserData.ax3);
    %currently can only update the following trial, since figure updates at
    %a time that isnt conducive to this?
    
    %This calculates latency to first lick
    if sum(scQtUserData.licks)>0
        if scQtUserData.soundDeliver(scQtUserData.trial) == 0 && scQtUserData.toneRule == 1 ||...
                scQtUserData.soundDeliver(scQtUserData.trial) == 1 && scQtUserData.toneRule == 0 ;
            scQtUserData.lickLat1(scQtUserData.lickLat1Counter,1) = scQtUserData.trial;
            scQtUserData.lickLat1(scQtUserData.lickLat1Counter,2) = scQtUserData.licks...
            (find(scQtUserData.licks>scQtUserData.soundOn(scQtUserData.trial),1,'first'))...
            -scQtUserData.soundOn(scQtUserData.trial);
            scQtUserData.lickLat1Counter = scQtUserData.lickLat1Counter + 1;
        elseif scQtUserData.soundDeliver(scQtUserData.trial) == 0 && scQtUserData.toneRule == 0 ||...
                scQtUserData.soundDeliver(scQtUserData.trial) == 1 && scQtUserData.toneRule == 1 ;
            scQtUserData.lickLat2(scQtUserData.lickLat2Counter,1) = scQtUserData.trial;
            scQtUserData.lickLat2(scQtUserData.lickLat2Counter,2) = scQtUserData.licks...
            (find(scQtUserData.licks>scQtUserData.soundOn(scQtUserData.trial),1,'first'))...
            -scQtUserData.soundOn(scQtUserData.trial);
            scQtUserData.lickLat2Counter = scQtUserData.lickLat2Counter + 1;
        end
        
    end
    
    if sum(scQtUserData.licks)>0
       if scQtUserData.soundDeliver(scQtUserData.trial) == 0 && scQtUserData.toneRule == 1 || scQtUserData.soundDeliver(scQtUserData.trial) == 1 && scQtUserData.toneRule == 0 ;
           scQtUserData.lickNum1(scQtUserData.lickNum1Counter,1) = scQtUserData.trial;
           scQtUserData.lickNum1(scQtUserData.lickNum1Counter,2) = size(find(...
            scQtUserData.licks>scQtUserData.soundOn(scQtUserData.trial) & ...
            scQtUserData.licks<scQtUserData.rewardTime(scQtUserData.trial)),1);
        scQtUserData.lickNum1Counter = scQtUserData.lickNum1Counter + 1;
       elseif scQtUserData.soundDeliver(scQtUserData.trial) == 0 && scQtUserData.toneRule == 0 || scQtUserData.soundDeliver(scQtUserData.trial) == 1 && scQtUserData.toneRule == 1 ;
           scQtUserData.lickNum2(scQtUserData.lickNum2Counter,1) = scQtUserData.trial;
           scQtUserData.lickNum2(scQtUserData.lickNum2Counter,2) = size(find(...
            scQtUserData.licks>scQtUserData.soundOn(scQtUserData.trial) & ...
            scQtUserData.licks<scQtUserData.rewardTime(scQtUserData.trial)),1);
        scQtUserData.lickNum2Counter = scQtUserData.lickNum2Counter + 1;
       end
    end
    
    %this prepares lick times for plotting as histogram
    if sum(scQtUserData.licks)>0
        scQtUserData.licks = round((scQtUserData.licks - scQtUserData.soundOn(scQtUserData.trial))/100)+20;
        scQtUserData.licks(scQtUserData.licks > 80) = [];
        scQtUserData.licks(scQtUserData.licks < 1) = [];
        if scQtUserData.soundDeliver(scQtUserData.trial) == 0 && scQtUserData.toneRule == 1 ...
                || scQtUserData.soundDeliver(scQtUserData.trial) == 1 && scQtUserData.toneRule == 0 ;
            scQtUserData.lickHist([scQtUserData.licks],1)=scQtUserData.lickHist([scQtUserData.licks],1)+1;
        elseif scQtUserData.soundDeliver(scQtUserData.trial) == 0 && scQtUserData.toneRule == 0 ...
                || scQtUserData.soundDeliver(scQtUserData.trial) == 1 && scQtUserData.toneRule == 1 ;
            scQtUserData.lickHist([scQtUserData.licks],2)=scQtUserData.lickHist([scQtUserData.licks],2)+1;
        end
    end
    
    
    title(scQtUserData.ax1,[...
        'Total trials: ', num2str(scQtUserData.trial),...
        ' Black = No Reward, Red = Reward']...
        ), ...
        
    %plot lick histogram
    plot(scQtUserData.lickAxes,scQtUserData.lickHist(:,1),'color','k','linewidth',2,...
        'parent',scQtUserData.ax1);
    plot(scQtUserData.lickAxes,scQtUserData.lickHist(:,2),'color','r','linewidth',2,...
        'parent',scQtUserData.ax1);
    
    %plot latency
    if sum(sum(scQtUserData.lickLat1)) == 0
        plot(scQtUserData.lickLat2(1:scQtUserData.lickLat2Counter-1,1),...
            scQtUserData.lickLat2(1:scQtUserData.lickLat2Counter-1,2),...
            'color','r','linestyle','none','marker','.','parent',scQtUserData.ax2)
    elseif sum(sum(scQtUserData.lickLat2)) == 0
        plot(scQtUserData.lickLat1(1:scQtUserData.lickLat1Counter-1,1),...
            scQtUserData.lickLat1(1:scQtUserData.lickLat1Counter-1,2),...
            'color','k','linestyle','none','marker','.','parent',scQtUserData.ax2)
    else
        plot(scQtUserData.lickLat1(1:scQtUserData.lickLat1Counter-1,1),...
            scQtUserData.lickLat1(1:scQtUserData.lickLat1Counter-1,2),...
            'color','k','linestyle','none','marker','.','parent',scQtUserData.ax2)
        plot(scQtUserData.lickLat2(1:scQtUserData.lickLat2Counter-1,1),...
            scQtUserData.lickLat2(1:scQtUserData.lickLat2Counter-1,2),...
            'color','r','linestyle','none','marker','.','parent',scQtUserData.ax2)
    end

    
%     %plot anticipatory licks
    if sum(sum(scQtUserData.lickNum1)) == 0
        plot(scQtUserData.lickNum2(1:scQtUserData.lickNum2Counter-1,1),...
            scQtUserData.lickNum2(1:scQtUserData.lickNum2Counter-1,2),...
            'color','r','linestyle','none','marker','.','parent',scQtUserData.ax3)
    elseif sum(sum(scQtUserData.lickNum2)) == 0
        plot(scQtUserData.lickNum1(1:scQtUserData.lickNum1Counter-1,1),...
            scQtUserData.lickNum1(1:scQtUserData.lickNum1Counter-1,2),...
            'color','k','linestyle','none','marker','.','parent',scQtUserData.ax3)
    else
        plot(scQtUserData.lickNum1(1:scQtUserData.lickNum1Counter-1,1),...
            scQtUserData.lickNum1(1:scQtUserData.lickNum1Counter-1,2),...
            'color','k','linestyle','none','marker','.','parent',scQtUserData.ax3)
        plot(scQtUserData.lickNum2(1:scQtUserData.lickNum2Counter-1,1),...
            scQtUserData.lickNum2(1:scQtUserData.lickNum2Counter-1,2),...
            'color','r','linestyle','none','marker','.','parent',scQtUserData.ax3)
    end
        
    scQtUserData.licks = zeros(1000,1);
    scQtUserData.lickCounter = 1;

    set(scQtUserData.ax1,'ygrid','on');
    set(scQtUserData.ax2,'ygrid','on');
    set(scQtUserData.ax3,'ygrid','on');
end

if (~isempty(strfind(newLine,'StartSession')))
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.ITI(scQtUserData.trial))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.soundRewDel(scQtUserData.trial))]);
    sendScQtControlMessage(['sound = ',num2str(scQtUserData.soundDeliver(scQtUserData.trial))]);
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'TriggerMatlab'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.ITI(scQtUserData.trial))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.soundRewDel(scQtUserData.trial))]);
    sendScQtControlMessage(['sound = ',num2str(scQtUserData.soundDeliver(scQtUserData.trial))]);
    sendScQtControlMessage('trigger(1)');
end

if (~isempty(strfind(newLine,'SoundOn'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.soundOn(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1));
end

if (~isempty(strfind(newLine,'Reward Delivered'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.rewardTime(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1));
end

if (~isempty(strfind(newLine,'Unrewarded'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.rewardTime(scQtUserData.trial) = str2num(newLine(1:spaceFinder(1)-1));
end

end