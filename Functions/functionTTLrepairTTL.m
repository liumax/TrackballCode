%this function is based on having a set number of target outputs, and
%trying to fix the actual TTLs

function [realTTLs] = functionTTLrepairTTL(expectedTTLs,realTTLs,threshDiffVal,maxDelayRatio,corrSlide);
disp('BEGINNING CHECK OF TTLS WITH EXPECTATIONS')
% threshDiffVal = 6;
% corrSlide = 100;
% maxDelayRatio = 1.2;

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
    bigFindOnset = find(onsetPhotDiff>max(expectedDiff)*maxDelayRatio);
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
        [xcf,lags,bounds]  = crosscorr(onsetPhotDiff,expectedDiff,corrSlide);
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
            %find longer thing, extend the other one to match
            if length(onsetPhotDiff) > length(expectedDiff)
                disp('RECEIVED TTLS GREATER THAN EXPECTED')
            elseif length(onsetPhotDiff) < length(expectedDiff)
                disp('RECEIVED TTLS LESS THAN EXPECTED')
                error('CANNOT COMPUTE')
            elseif length(onsetPhotDiff) == length(expectedDiff)
                disp('RECEIVED TTLS EQUAL TO EXPECTED')
            end
%             diffVal = abs(reshape(onsetPhotDiff,1,[])-reshape(expectedDiff,1,[]));
            whileTrig = 0;
            crawlInd = 1;
            disp('BEGIN DATA CRAWL')
            while whileTrig == 0;
%                 disp(crawlInd)
                diffVal = onsetPhotDiff(crawlInd) - expectedDiff(crawlInd);
                if diffVal > threshDiffVal
                    realTTLs(crawlInd) = [];
                    onsetPhotDiff = diff(realTTLs);
                    if length(realTTLs) == length(expectedTTLs)
                        disp('REMOVING POINT WITH CRAWLER, PROBLEM FIXED')
                        crawlStore = crawlInd;
                        break
                    else
                        disp('PROBLEM NOT SOLVED WITH DATA CRAWLER FULLY, CONTINUING...')
                        crawlInd = crawlInd + 1;
                    end
                elseif diffVal < -threshDiffVal
                    realTTLs(crawlInd+1) = [];
                    onsetPhotDiff = diff(realTTLs);
                    if length(realTTLs) == length(expectedTTLs)
                        disp('REMOVING POINT WITH CRAWLER, PROBLEM FIXED')
                        crawlStore = crawlInd;
                        break
                    else
                        disp('PROBLEM NOT SOLVED WITH DATA CRAWLER FULLY, CONTINUING...')
                        crawlInd = crawlInd + 1;
                    end
                else
                    crawlInd = crawlInd + 1;
                end
                if crawlInd >= length(expectedDiff)
                    disp('END OF CRAWL')
                    break
                end
            end
        end
    end
end


end