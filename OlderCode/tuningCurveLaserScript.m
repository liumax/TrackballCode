%This is the Matlab Script

global scQtUserData;
          
% UI prompt:
prompt = {'Mouse ID:',...  
    'Left/Right:',...       
    'Recording Depth',...
    'Laser Pulse Number:',...           
    'Laser Pulse Width:',...           
    'Laser Pulse ITI:',...        
    'Notes:'}; %the bracket is to end the prompt     
dlg_title = 'tuningCurveLaser:';
num_lines=1;
def={'','','','20','10','40',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
pause(2); % need to pause for microcontroller or things break!

% Connect to Microcontroller and store MetaData:
% setComVal; % Call script to set com port string for this computer
% sHandle = scConnect(comValStr,@pokeE_Callback);
% pause(1);

i=1;
t = clock;
rand('seed',sum(round(clock)));
scQtUserData.taskID = 'tuningCurveLaser';
scQtUserData.mouseID = answer{i};i=i+1;
scQtUserData.side = answer{i};i=i+1;
scQtUserData.depth = str2num(answer{i});i=i+1;
scQtUserData.laserNum = str2num(answer{i});i=i+1;
scQtUserData.laserDur = str2num(answer{i});i=i+1;
scQtUserData.laserITI = str2num(answer{i});i=i+1;
scQtUserData.notes = answer{i};i=i+1;

%Date/Time
scQtUserData.date = date;
scQtUserData.time = strcat(num2str(t(4)),':',num2str(t(5)));

pause(0.2);

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''side:', (scQtUserData.side),''')']);
sendScQtControlMessage(['disp(''depth:', num2str(scQtUserData.depth),''')']);
sendScQtControlMessage(['disp(''laserNum:', num2str(scQtUserData.laserNum),''')']);
sendScQtControlMessage(['disp(''laserDur:', num2str(scQtUserData.laserDur),''')']);
sendScQtControlMessage(['disp(''laserITI:', num2str(scQtUserData.laserITI),''')']);
sendScQtControlMessage(['disp(''taskID:', scQtUserData.taskID,''')']);
sendScQtControlMessage(['disp(''date:', scQtUserData.date,''')']);
sendScQtControlMessage(['disp(''time:', scQtUserData.time,''')']);
sendScQtControlMessage(['disp(''notes:', scQtUserData.notes,''')']);

pause(1) %Need to put all my timings in before this stuff


sendScQtControlMessage(['laserNum=',num2str(scQtUserData.laserNum)]);
sendScQtControlMessage(['laserDur=',num2str(scQtUserData.laserDur)]);
sendScQtControlMessage(['laserITI=',num2str(scQtUserData.laserITI)]);
