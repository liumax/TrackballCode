global scQtUserData

pause(1);

% sendScQtControlMessage('portout[4] = 0');
% sendScQtControlMessage('portout[5] = 0');

% UI prompt:
prompt = {'Side (r/l/b):',... 1
    'Num Iter:',...  2
    'ITI (ms)',... 3
    'OpenDur (ms):'}; %          4
dlg_title = 'pokeF:';

answer = inputdlg(prompt,dlg_title);
t = clock;

% Connect to Microcontroller and store MetaData:
scQtUserData.side   = answer{1};
scQtUserData.numIter = str2double(answer{2});
scQtUserData.iti   = str2double(answer{3});
scQtUserData.openDur    = str2double(answer{4});

pause(0.2);
sendScQtControlMessage(['openDur = ', num2str(scQtUserData.openDur)]);

for trialNum = 1:scQtUserData.numIter
    sendScQtControlMessage(['disp(''Trial ', num2str(trialNum),''')']);
    if strcmp(scQtUserData.side,'l')
        sendScQtControlMessage('trigger(1)');
    elseif strcmp(scQtUserData.side,'r')
        sendScQtControlMessage('trigger(2)');
    elseif strcmp(scQtUserData.side,'b')
        sendScQtControlMessage('trigger(3)');
    end
    pause((scQtUserData.iti - scQtUserData.openDur)./1000);
end