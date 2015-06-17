%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

int upA = 0
int velCount = 0

int intWindow1 = 500

callback portin[7] up 
    upA = upA + 1
    if upA == 3 do
        velCount = velCount + 1
        disp(velCount)
        upA = 0
        do in intWindow1
            velCount = velCount - 1
            disp(velCount)
        end
    end
end;
