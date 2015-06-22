%This is the Matlab Script

global scQtUserData;

scQtUserData.mouseID = 'HelloKitty';
scQtUserData.tester = 100;
scQtUserData.minRew = 200;

pause(1)

sendScQtControlMessage(['disp(''Mouse ID:', scQtUserData.mouseID,''')']);
sendScQtControlMessage(['disp(''tester ID:', num2str(scQtUserData.tester),''')']);
sendScQtControlMessage(['disp(''minRew:', num2str(scQtUserData.minRew),''')']);
