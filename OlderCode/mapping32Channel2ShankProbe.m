%the array for displaying how things are should be 20 x 4. This is divided
%into a top and bottom set of 10 x 4.

myShankDesig = 1:1:32;

neuroTechShankDesig = [3
16
5
12
9
8
11
6
7
2
1
4
13
14
10
15
29
17
28
21
30
22
32
25
31
27
26
24
20
19
23
18];


neuroTechDesig = [...
    11,0,0,32;...
    9,0,0,30;...
    7,0,0,31;...
    5,0,0,28;...
    3,1,26,29;...
    2,4,24,27;...
    6,13,20,25;...
    8,14,19,22;...
    10,15,18,23;...
    12,16,17,21];

neuroNexusDesig = [...
    1,0,0,32;...
    2,0,0,31;...
    3,0,0,30;...
    4,0,0,29;...
    5,16,17,28;...
    6,15,18,27;...
    7,14,19,26;...
    8,13,20,25;...
    9,12,21,24;...
    10,11,22,23];
%This is looking down at omnetics connector from the top. The top omnetics
%connector is the one on the front of the PCB. 
neuroNexusOmnetics = [...
    23,25,27,29,31,19,17,21,11,15,13,1,3,5,7,9;...
    24,26,28,30,32,20,18,22,12,16,14,2,4,6,8,10];
%This is the trodes designations looking at the pins. I realize now that
%this is actually the mirror image of what I should want. Will repair that
%below. 
trodesOmneticsSingle = [...
    8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23;
    7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24];
trodesOmneticsSingle = flipud(fliplr(trodesOmneticsSingle));

trodesOmnetics = trodesOmneticsSingle;

%generate master sheet that will hold everything
master = zeros(length(myShankDesig),6);
master(:,1) = myShankDesig;
master(:,2) = neuroTechShankDesig;
%now I want to find the appropriate neuronexus numbers for the same
%positions
for i = 1:32
    posFinder = find(neuroTechDesig==master(i,2));
    master(i,3) = neuroNexusDesig(posFinder);
end
%now convert this to the omnetics connector to get trodes hardware
%designations
for i = 1:32
    posFinder = find(neuroNexusOmnetics == master(i,3));
    master(i,4) = trodesOmnetics(posFinder);
end



