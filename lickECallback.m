%This is the Matlab Callback

function lickACallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    lickEScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,'Trial = 400')))
    scQtUserData.tripSwitch = 1;
end


if scQtUserData.trial>1 && (~isempty(strfind(newLine,'SoundOff')))
    scQtUserData.master(scQtUserData.trial,11) = (scQtUserData.master(scQtUserData.trial,10)...
        -scQtUserData.master(scQtUserData.trial-1,10))/1000;
    sendScQtControlMessage(['disp(''ITI = ',num2str(scQtUserData.master(scQtUserData.trial,11)),''')']);
    
    if ~isfield(scQtUserData,'updateFig')
        disp('resetting updateFig');
        scQtUserData.updateFig = -1;
    end
    if ~ishandle(scQtUserData.updateFig)
        scQtUserData.updateFig = figure('color','w');
        scQtUserData.corrAx = subplot(2,1,1,'parent',scQtUserData.updateFig);
        scQtUserData.durAx = subplot(2,1,2,'parent',scQtUserData.updateFig);
        hold(scQtUserData.corrAx,'on');
        hold(scQtUserData.durAx,'on');
        ylabel(scQtUserData.corrAx,'Licks');
        ylabel(scQtUserData.durAx,'Reward Dur (ms)');
    end

    cla(scQtUserData.corrAx);
    
    title(scQtUserData.corrAx,[...
        'Total trials ', num2str(scQtUserData.trial),...
        ' Black = prelick',' Blue = Anticipatory', ' Green = post', 'Red = other']), ...
    % Plot pre-licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,6),'linewidth',2,...
        'color','k','parent',scQtUserData.corrAx);
    % Plot anticipatory licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,7),'linewidth',2,...
        'color','b','parent',scQtUserData.corrAx);
    %Plot post-delivery licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,8),'linewidth',2,...
        'color','g','parent',scQtUserData.corrAx);
    %Plot all other licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,9),'linewidth',2,...
        'color','r','parent',scQtUserData.corrAx);
    %Plot all other licks
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,11),'linewidth',3,...
        'color','c','parent',scQtUserData.corrAx);
    % Plot reward size
    plot(1:scQtUserData.trial,scQtUserData.master(1:scQtUserData.trial,1),'color','k','linewidth',3,...
        'parent',scQtUserData.durAx);

    set(scQtUserData.corrAx,'ylim',[0 30],'ytick',0:1:30,'ygrid','on');
    set(scQtUserData.durAx,'ygrid','on');
end

if (~isempty(strfind(newLine,'SoundOff'))) && scQtUserData.tripSwitch == 0;
    scQtUserData.trial = scQtUserData.trial + 1;
    sendScQtControlMessage(['disp(''Trial = ',num2str(scQtUserData.trial),''')']);
    sendScQtControlMessage(['itiDur = ',num2str(scQtUserData.master(scQtUserData.trial,2))]); 
    sendScQtControlMessage(['soundRewDel = ',num2str(scQtUserData.master(scQtUserData.trial,3))]);
    sendScQtControlMessage(['rewLength = ',num2str(scQtUserData.master(scQtUserData.trial,1))]);
    spaceFinder = find(newLine == ' ');
    scQtUserData.master(scQtUserData.trial,10) = str2num(newLine(1:spaceFinder(1)-1));
%     sendScQtControlMessage(['disp(''SoundTime = ',num2str(scQtUserData.master(...
%         scQtUserData.trial,10)),''')']);
    sendScQtControlMessage('trigger(1)');
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


end