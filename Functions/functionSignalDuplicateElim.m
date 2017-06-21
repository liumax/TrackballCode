%This is a function for the elimination of duplicate values, which commonly
%occur in MBED outputs, since they save every triggered event.

%Inputs: dataInput should be data organized in columns. 
%columnInd should be the index that I will be looking for duplicates in.

%Outputs: revData: revised data with eliminated duplicates

function [revData] = functionSignalDuplicateElim(dataInput,columnInd);

whileTrig = 1;
whileInd = 2;


while whileTrig == 1;
    if whileInd > length(dataInput)
        break
    end
    prevVal = dataInput(whileInd-1,columnInd);
    currVal = dataInput(whileInd,columnInd);
    if prevVal == currVal
%         disp('Duplicate Detected!')
        dataInput(whileInd,:) = [];
    else
        whileInd = whileInd + 1;
    end
end

revData = dataInput;

end