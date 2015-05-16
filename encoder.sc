%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

int upA = 0
int upB = 0
int downA = 0
int downB = 0

int intWindow1 = 500

callback portin[7] up 
    upA = upA + 1
    do in intWindow1
        upA = upA - 1
    end
    disp(upA)
end;
