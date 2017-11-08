%this function is based on having a set number of target outputs, and
%trying to fix the actual TTLs

function [realTTLs] = functionTTLrepairTTL(expectedTTLs,realTTLs);
disp('BEGINNING CHECK OF TTLS WITH EXPECTATIONS')

%first, calculate differences
onsetPhotDiff = diff(realTTLs);
expectedDiff = diff(expectedTTLs);

if length(realTTLs) == length(expectedTTLs);
    disp('Actual Pulses Equal to Planned Number of Pulses, Proceeding')
else
    %calculate actual ITIs
    
    disp('Actual Pulses NOT Equal To Planned Pulses')
    disp(strcat('Actual Pulses:',num2str(length(realTTLs))))
    disp(strcat('Expected Pulses:',num2str(length(expectedTTLs))))
    %test of there are delays in onset phot that are huge
    bigFindOnset = find(onsetPhotDiff>max(expectedDiff)*1200);
    if bigFindOnset
        disp('FOUND EXTRA LONG PULSE, DELETING')
        realTTLs(bigFindOnset) = [];
        onsetPhotDiff = diff(realTTLs);
        %check if this fixes length issue
        if length(realTTLs) == length(expectedTTLs)
            disp('LENGTH ERROR FIXED')
        else
            error('PROBLEM NOT FIXED')
        end
    else
        disp('NO BIG ITIS FOUND...LOOKING FOR MISMATCH WITH CROSSCORR')
        [xcf,lags,bounds]  = crosscorr(onsetPhotDiff/1000,expectedDiff,300);
        [xcMax maxInd] = max(xcf);
        xcLag = lags(maxInd);
        disp('CorrLag')
        disp(xcLag)
        disp('MaxCorr')
        disp(xcMax)
        if xcMax > 0.9
            if xcLag < 0
                disp('DELETING EXCESS PULSES')
                realTTLs(1:-xcLag) = [];
                onsetPhotDiff = diff(realTTLs);
                %check if lengths are right
                if length(realTTLs) == length(expectedTTLs)
                    disp('FIXED THE PROBLEM')
                else
                    error('PROBLEM NOT FIXED')
                end
            elseif xcLag > 0
                error('ALIGNMENT IS IN THE POSITIVE DIRECTION...ERROR')
            elseif xcLag == 0
                error('TTLS ALIGN WITH DELAYS')
            end
        else
            %time to crawl
            diffVal = abs(onsetPhotDiff-expectedDiff);
            whileTrig = 0;
            crawlInd = 1;
            while whileTrig == 0;
                if diffVal(crawlInd) > 0.6
                    realTTLs(crawlInd) = [];
                    onsetPhotDiff = diff(realTTLs);
                    if length(realTTLs) == length(expectedTTLs)
                        disp('PROBLEM SOLVED WITH DATA CRAWLER')
                        crawlStore = crawlInd;
                        break
                    else
                        error('PROBLEM NOT SOLVED WITH DATA CRAWLER, KILLING PROGRAM')
                    end
                else
                    crawlInd = crawlInd + 1;
                end
                if crawlInd >= length(diffVal)
                    break
                end
            end
        end
    end
end


end