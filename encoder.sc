%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

callback portin[7] up % lickometer activated
    disp('InputA')
end;

callback portin[8] up % lickometer activated
    disp('InputB')
end;
