%This is meant to be a basic function for finding strings in a cell array.
%The goal of this function is to input a cell array and the string you are
%interested in, and get an output of the index of the string. 

function [findString] = functionCellStringFind(cellArray,string);

%find the string
findString = strfind(cellArray,string);
%convert to something useable with cellfun
findString = cellfun('isempty',findString);
%now find if a point exists
findString = find(findString == 0);

if isempty(findString)
    disp('Target String Not Found')
end

end

