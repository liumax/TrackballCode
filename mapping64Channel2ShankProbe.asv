%the array for displaying how things are should be 20 x 4. This is divided
%into a top and bottom set of 10 x 4.

myShankDesig = 1:1:64;

neuroTechShankDesig = [21
23
24
30
29
16
18
20
27
19
17
25
26
32
28
22
1
3
5
7
9
11
13
15
31
14
12
10
8
6
4
2
64
62
60
58
56
54
52
50
34
51
53
55
57
59
61
63
44
42
41
35
36
49
47
45
38
46
48
40
39
33
37
43];


neuroTechDesig = [...
    1,0,0,64;...
    2,0,0,63;...
    3,0,0,62;...
    4,0,0,61;...
    5,6,59,60;...
    7,8,57,58;...
    9,10,55,56;...
    11,12,53,54;...
    13,14,51,52;...
    15,16,49,50;...
    17,0,0,48;...
    18,0,0,47;...
    19,0,0,46;...
    20,0,0,45;...
    21,22,43,44;...
    23,28,37,42;...
    24,32,33,41;...
    25,29,36,40;...
    26,30,35,39;...
    27,31,34,38];

neuroNexusDesig = [...
    24,0,0,41;...
    27,0,0,38;...
    22,0,0,43;...
    28,0,0,37;...
    29,26,39,36;...
    20,25,40,45;...
    18,23,42,47;...
    16,21,44,49;...
    14,19,46,51;...
    12,17,48,53;...
    10,0,0,55;...
    15,0,0,50;...
    8,0,0,57;...
    13,0,0,52;...
    6,11,54,59;...
    4,5,60,61;...
    9,32,33,56;...
    2,3,62,63;...
    7,30,35,58;...
    1,31,34,64];
%This is looking down at omnetics connector from the top. The top omnetics
%connector is the one on the front of the PCB. 
neuroNexusOmnetics = [...
    34,35,62,33,60,54,57,55,10,8,11,5,32,3,30,31;...
    64,58,63,56,61,59,52,50,15,13,6,4,9,2,7,1;...
    53,51,49,47,45,36,37,38,27,28,29,20,18,16,14,12;...
    48,46,44,42,40,39,43,41,24,22,26,25,23,21,19,17];
%This is the trodes designations looking at the pins. I realize now that
%this is actually the mirror image of what I should want. Will repair that
%below. 
trodesOmneticsSingle = [...
    8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23;
    7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24];
% trodesOmneticsSingle = fliplr(trodesOmneticsSingle); %actually realized I
% had the right orientation all along. XD
%generate the inverse for the flipped headstage. add 
trodesOmneticsFlip = rot90(rot90(trodesOmneticsSingle))+32;
%combine the two headstages
trodesOmnetics = [trodesOmneticsSingle;trodesOmneticsFlip];

%generate master sheet that will hold everything
master = zeros(length(myShankDesig),6);
master(:,1) = myShankDesig;
master(:,2) = neuroTechShankDesig;
%now I want to find the appropriate neuronexus numbers for the same
%positions
for i = 1:64
    posFinder = find(neuroTechDesig==master(i,2));
    master(i,3) = neuroNexusDesig(posFinder);
end
%now convert this to the omnetics connector to get trodes hardware
%designations
for i = 1:64
    posFinder = find(neuroNexusOmnetics == master(i,3));
    master(i,4) = trodesOmnetics(posFinder);
end



