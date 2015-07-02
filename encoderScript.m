%This is the Matlab Script

global scQtUserData;


scQtUserData.velocity = zeros(1000000,2);
scQtUserData.velCounter = 1; %This is to maintain correct position in velocity
scQtUserData.plotCounter = 0; %This is to trigger plotting