fname = '2msDiff100xWithPulsing2msMCU4msAuditoryTTL';
%inputs portstates
[portStates] = maxTrialVariableNoTask(fname);
%finds portstate upswings
stateInputs = diff(portStates.inStates(:,2));
finder = find(stateInputs == 1)+1;
finder = [1;finder];
%finds the times of these upswings
times = portStates.tStamps(finder);
timeDiff= diff(times);

MCU = find(diff(portStates.outStates(:,1))==1);
Laser = find(diff(portStates.outStates(:,3))==1);

stateOutputs = diff(portStates.outStates(:,));

diffSmall= find(timeDiff<12);
diffLarge = find(timeDiff>12);

plot(timeDiff)