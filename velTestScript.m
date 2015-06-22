%This is the Matlab Script

global scQtUserData;


scQtUserData.velocity = zeros(1,2);
scQtUserData.velCounter = 1; %This is to maintain correct position in velocity
scQtUserData.plotCounter = 0; %This is to trigger plotting
scQtUserData.firstPoint = 1; %This is a placeholder for the first point to plot