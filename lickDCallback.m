%This is the Matlab Callback

function lickACallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    lickDScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,'Trial = 400')))
    scQtUserData.tripSwitch = 1;
end

if scQtUserData.trial>=1 && (~isempty(strfind(newLine,'TriggerMatlab')))
    if scQtUserData.trial>1
        scQtUserData.master(scQtUserData.trial,11) = (scQtUserData.master(scQtUserData.trial,10)...
            -scQtUserData.master(scQtUserData.trial-1,10))/1000; %This is to log sound off times and ITIs
%     sendScQtControlMessage(['disp(''ITI = ',num2str(scQtUserData.master(scQtUserData.trial,11)),''')']);
    end
    
    if ~isfield(scQtUserData,'updateFig') %This code is just to make sure updateFig has a value.
        disp('resetting updateFig');
        scQtUserData.updateFig = -1;
    end
    
    if ~ishandle(scQtUserData.updateFig) %This is to set the basis for all the plots!
        scQtUserData.updateFig = figure('color','w');
        scQtUserData.lickAx = subplot(3,1,1,'parent',scQtUserData.updateFig);
        scQtUserData.durAx = subplot(3,1,2,'parent',scQtUserData.updateFig);
        scQtUserData.velAx = subplot(3,1,3,'parent',scQtUserData.updateFig);
        hold(scQtUserData.lickAx,'on');
        hold(scQtUserData.durAx,'on');
        hold(scQtUserData.velAx,'on');
        ylabel(scQtUserData.lickAx,'Licks');
        ylabel(scQtUserData.durAx,'Reward Dur (ms)');
        ylabel(scQtUserData.velAx,'Speed (cm/s)');
    end

    cla(scQtUserData.lickAx); %This clears lickAx. This is because there are constantly adjustments to lickAx
    
    %currently can only update the following trial, since figure updates at
    %a time that isnt conducive to this?
%     scQtUserData.velStart(scQtUserData.trial) = find(scQtUserData.velocity(:,1)<scQtUserData.master(scQtUserData.trial,10)-5000);
%     scQtUserData.velEnd(scQtUserData.trial) = find(scQtUserData.velocity(:,1)>scQtUserData.master(scQtUserData.trial,10)+2000);
    
    velFinder = find(scQtUserData.velocity(:,1)>scQtUserData.master(scQtUserData.trial,10)-5000 & ...
        scQtUserData.velocity(:,1)<scQtUserData.master(scQtUserData.trial,10)+2000);
    if ~isempty(velFinder)
        velStart = velFinder(1)
        velEnd = velFinder(end)
        velEnd-velStart
    end
    
    
    title(scQtUserData.lickAx,[...
        'Total trials ', num2str(scQtUserData.trial),...
        ' Black = prelick',' Blue = Anticipatory', ' Green = post', 'Red = other']), ...
    % Plot pre-licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,6),'linewidth',2,...
        'color','k','parent',scQtUserData.lickAx);
    % Plot anticipatory licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,7),'linewidth',2,...
        'color','b','parent',scQtUserData.lickAx);
    %Plot post-delivery licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,8),'linewidth',2,...
        'color','g','parent',scQtUserData.lickAx);
    %Plot all other licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,9),'linewidth',2,...
        'color','r','parent',scQtUserData.lickAx);
    %Plot all other licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,11),'linewidth',3,...
        'color','c','parent',scQtUserData.lickAx);
    % Plot reward size
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,1),'color','k','linewidth',3,...
        'parent',scQtUserData.durAx);
    %Plot out velocity during trial!
    size(scQtUserData.velocity(velStart,1)...
            :scQtUserData.velocity(velEnd,1))
    size(scQtUserData.velocity(velStart,2)...
            :scQtUserData.velocity(velEnd,2)) 
    
%     if ~isempty(velStart) && ~isempty(velEnd) 
%         plot((scQtUserData.velocity(velStart,1)...
%             :scQtUserData.velocity(velEnd,1)),...
%             scQtUserData.velocity(velStart,2):...
%             scQtUserData.velocity(velEnd,2),...
%             'color','k','linewidth',2,'parent',scQtUserData.velAx);
%     end

    set(scQtUserData.lickAx,'ylim',[0 30],'ytick',0:1:30,'ygrid','on');
    set(scQtUserData.durAx,'ygrid','on');
    set(scQtUserData.velAx,'ygrid','on');
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

if (~isempty(strfind(newLine,'SoundOff'))) && scQtUserData.tripSwitch == 0;
    spaceFinder = find(newLine == ' ');
    scQtUserData.master(scQtUserData.trial,10) = str2num(newLine(1:spaceFinder(1)-1));
end

if (~isempty(strfind(newLine,'trialState')))
    findSpaces=find(newLine == ' ');
    stateValue=str2double(newLine(findSpaces(end)+1:end));
    disp(stateValue)
    if stateValue == 1
        scQtUserData.master(scQtUserData.trial,6)=scQtUserData.master(scQtUserData.trial,6)+1;
%         disp '1'
        disp(scQtUserData.master(scQtUserData.trial,6))
    elseif stateValue == 2
        scQtUserData.master(scQtUserData.trial,7)=scQtUserData.master(scQtUserData.trial,7)+1;
%         disp '2'
        disp(scQtUserData.master(scQtUserData.trial,7))
    elseif stateValue == 3
        scQtUserData.master(scQtUserData.trial,8)=scQtUserData.master(scQtUserData.trial,8)+1;
%         disp '3'
        disp(scQtUserData.master(scQtUserData.trial,8))
    elseif stateValue == 4
        scQtUserData.master(scQtUserData.trial,9)=scQtUserData.master(scQtUserData.trial,9)+1;
%         disp '4'
        disp(scQtUserData.master(scQtUserData.trial,9))
    end
    disp 'State Value Updated'
end

if (~isempty(strfind(newLine,'upA')))
    findSpacer=find(newLine == ' ');
    scQtUserData.velocity(scQtUserData.velCounter,1) = str2double(newLine(1:(findSpacer(1)-1)));
    scQtUserData.velocity(scQtUserData.velCounter,2) = str2double(newLine(findSpacer(end)+1:end));
    x=scQtUserData.velocity(scQtUserData.velCounter,1)
    y=scQtUserData.velocity(scQtUserData.velCounter,2)
    scQtUserData.velCounter = scQtUserData.velCounter + 1;
end

end