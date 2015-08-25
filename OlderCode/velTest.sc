%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

int upA = 0
int intWindow = 500

callback portin[7] up 
    upA = upA + 1
    disp(upA)
    do in intWindow
        upA = upA - 1
        disp(upA)
    end
end;
