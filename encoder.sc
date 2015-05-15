%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

int upA = 0
int upB = 0
int downA = 0
int downB = 0

callback portin[7] up 
    upA = upA + 1
end;

callback portin[7] down 
    downA = downA + 1
end;

callback portin[8] up 
    upB = upB + 1
end;

callback portin[8] down 
    downB = downB + 1
end;

