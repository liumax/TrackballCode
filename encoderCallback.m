%This is the Matlab Callback

function encoderCallback(newLine)

global scQtHistory; %multipurpose place to store processed event history
global scQtControllerOutput; %the text output from the microcontroller 
global scQtCallBackHandle; %the handle to the function called for every new event
global scQtInitiated; %the callback function should set this to 1 once all user variables are set

global scQtUserData;

if (scQtInitiated == 0)
    encoderScript;
    scQtInitiated = 1;
    newLine = 'start next trial';
end

if (~isempty(strfind(newLine,'upA')))
    findSpacer=find(newLine == ' ');
    scQtUserData.velocity(scQtUserData.velCounter,1) = str2double(newLine(1:(findSpacer(1)-1)));
    scQtUserData.velocity(scQtUserData.velCounter,2) = str2double(newLine(findSpacer(end)+1:end));
    scQtUserData.velCounter = scQtUserData.velCounter + 1;
    scQtUserData.plotCounter = scQtUserData.plotCounter + 1;
end

if scQtUserData.plotCounter == 10
    scQtUserData.plotCounter = 0;
    if ~isfield(scQtUserData,'updateFig') %This code is just to make sure updateFig has a value.
        disp('resetting updateFig');
        scQtUserData.updateFig = -1;
    end
    
    if ~ishandle(scQtUserData.updateFig) %This is to set the basis for all the plots!
        scQtUserData.updateFig = figure('color','w');
        scQtUserData.velAx = subplot(2,1,1,'parent',scQtUserData.updateFig);
        hold(scQtUserData.velAx,'on');
        xlabel(scQtUserData.velAx,'Time (s)');
        ylabel(scQtUserData.velAx,'Speed (cm/s)');
    end
    
    cla(scQtUserData.velAx); %This clears lickAx. This is because there are constantly adjustments to lickAx

    %currently can only update the following trial, since figure updates at
    %a time that isnt conducive to this?
%     scQtUserData.velStart(scQtUserData.trial) = find(scQtUserData.velocity(:,1)<scQtUserData.master(scQtUserData.trial,10)-5000);
%     scQtUserData.velEnd(scQtUserData.trial) = find(scQtUserData.velocity(:,1)>scQtUserData.master(scQtUserData.trial,10)+2000);
    
    title(scQtUserData.velAx,'Mouse Speed (cm/s)'), ...
    plot((scQtUserData.velocity(:,1)-scQtUserData.velocity(1,1))/1000,scQtUserData.velocity(:,2)*2*0.1,'color','k','linewidth',2,...
    'parent',scQtUserData.velAx);

    set(scQtUserData.velAx,'ygrid','on');
end

end