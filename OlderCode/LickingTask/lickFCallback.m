%This is the Matlab Callback

function lickACallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    lickFScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,['Trial = ',num2str(scQtUserData.blocks*scQtUserData.blockSize)])))
    scQtUserData.tripSwitch = 1;
end

if scQtUserData.trial>=1 && (~isempty(strfind(newLine,'TriggerMatlab')))
    if scQtUserData.trial>1
        scQtUserData.master(scQtUserData.trial,11) = (scQtUserData.master(scQtUserData.trial,10)...
            -scQtUserData.master(scQtUserData.trial-1,10))/1000; %This is to log sound off times and ITIs
    end
    
    if ~isfield(scQtUserData,'updateFig') %This code is just to make sure updateFig has a value.
        disp('resetting updateFig');
        scQtUserData.updateFig = -1;
    end
    
    if ~ishandle(scQtUserData.updateFig) %This is to set the basis for all the plots!
        scQtUserData.updateFig = figure('color','w');
        scQtUserData.lickAx = subplot(3,1,1,'parent',scQtUserData.updateFig);
        scQtUserData.durAx = subplot(3,1,2,'parent',scQtUserData.updateFig);
        scQtUserData.histAx = subplot(3,1,3,'parent',scQtUserData.updateFig);
        hold(scQtUserData.lickAx,'on');
        hold(scQtUserData.durAx,'on');
        hold(scQtUserData.histAx,'on');
        ylabel(scQtUserData.lickAx,'Time (s)');
        ylabel(scQtUserData.durAx,'Reward Dur (ms)');
        ylabel(scQtUserData.histAx,'Licks');
    end

    cla(scQtUserData.lickAx); %This clears lickAx. This is because there are constantly adjustments to lickAx
    cla(scQtUserData.histAx);
    %currently can only update the following trial, since figure updates at
    %a time that isnt conducive to this?
    
    
    scQtUserData.licks = round((scQtUserData.licks - scQtUserData.master(scQtUserData.trial,10))/100)+20;
    scQtUserData.licks(scQtUserData.licks > 80) = [];
    scQtUserData.licks(scQtUserData.licks < 1) = [];
    
    if scQtUserData.master(scQtUserData.trial,1) == scQtUserData.minRew;
        scQtUserData.lickHist([scQtUserData.licks],1)=scQtUserData.lickHist([scQtUserData.licks],1)+1;
    elseif scQtUserData.master(scQtUserData.trial,1) == scQtUserData.maxRew;
        scQtUserData.lickHist([scQtUserData.licks],2)=scQtUserData.lickHist([scQtUserData.licks],2)+1;
    end
    
    
    title(scQtUserData.lickAx,[...
        'Total trials: ', num2str(scQtUserData.trial),...
        ' Completed trials: ', num2str(sum(scQtUserData.master(:,6))),...
        ' Incomplete Trials: ', num2str(sum(scQtUserData.master(:,7))),...
        ' Premature Licking: ', num2str(sum(scQtUserData.master(:,8))),...
        ' C = ITI']), ...

    %Plot ITI
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,11),'linewidth',3,...
        'color','c','parent',scQtUserData.lickAx);
    % Plot reward size
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,1),'color','k','linewidth',3,...
        'parent',scQtUserData.durAx);
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,12),'color','g','linestyle','none',...
        'marker','.','parent',scQtUserData.durAx)
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,13),'color','r','linestyle','none',...
        'marker','.','parent',scQtUserData.durAx)
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,14),'color','b','linestyle','none',...
        'marker','o','parent',scQtUserData.durAx)
    plot(scQtUserData.lickAxes,scQtUserData.lickHist(:,1),'color','k','linewidth',2,...
        'parent',scQtUserData.histAx);
    plot(scQtUserData.lickAxes,scQtUserData.lickHist(:,2),'color','r','linewidth',2,...
        'parent',scQtUserData.histAx);
%     
    
    scQtUserData.licks = zeros(1000,1);
    scQtUserData.lickCounter = 1;

    set(scQtUserData.lickAx,'ygrid','on');
    set(scQtUserData.durAx,'ygrid','on');
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

if ~isempty(strfind(newLine,'Lick Detected')) && scQtUserData.tripSwitch == 0;
    spaceFinder= find(newLine == ' ');
    scQtUserData.licks(scQtUserData.lickCounter,1) = str2num(newLine(1:spaceFinder(1)-1));
    scQtUserData.lickCounter= scQtUserData.lickCounter+1;
end

if (~isempty(strfind(newLine,'SoundOn'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.master(scQtUserData.trial,10) = str2num(newLine(1:spaceFinder(1)-1));
end

if (~isempty(strfind(newLine,'Reward Delivered'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.master(scQtUserData.trial,6) = 1;
    scQtUserData.master(scQtUserData.trial,12) = scQtUserData.master(scQtUserData.trial,1);
end

if (~isempty(strfind(newLine,'No Reward'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.master(scQtUserData.trial,7) = 1;
    scQtUserData.master(scQtUserData.trial,13) = scQtUserData.master(scQtUserData.trial,1);
end

if (~isempty(strfind(newLine,'Bad Licks'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.master(scQtUserData.trial,8) = 1;
    scQtUserData.master(scQtUserData.trial,14) = scQtUserData.master(scQtUserData.trial,1);
end


end