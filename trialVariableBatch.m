clear all

%Note: ls only outputs a matrix of m x n in WINDOWS NOT IN MACINTOSH
filelist=ls('*.txt');
[m n]=size(filelist);


for i=1:m
    fname=filelist(i,:);
    [trialStates, portStates, trialParams] = maxTrialVariables(fname);
    [cleanLicks] = maxTrialVariableAnalysis(trialStates,trialParams,portStates);
    holder{i} = cleanLicks;
end




   