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

neuraLynxDesig = [...
    0,NaN,NaN,63;...
    1,NaN,NaN,62;...
    2,NaN,NaN,61;...
    3,NaN,NaN,60;...
    4,5,58,59;...
    6,7,56,57;...
    8,9,54,55;...
    10,11,52,53;...
    12,13,50,51;...
    14,15,48,49;...
    16,NaN,NaN,47;...
    17,NaN,NaN,46;...
    18,NaN,NaN,45;...
    19,NaN,NaN,44;...
    20,21,42,43;...
    22,27,36,41;...
    23,31,32,40;...
    24,28,35,39;...
    25,29,34,38;...
    26,30,33,37];
%This is looking down at omnetics connector from the top. Here, for
%convenience I have continued the old system of having the connectors
%stacked on top of each other, even though they are actually side by side.
%Here, up in the matrix indicates the surface away from the actual PCB. 
neuralynxOmnetics = [...
    31:-1:16;...
    15:-1:0;...
    63:-1:48;...
    47:-1:32];
%This is the trodes designations looking at the pins. I realize now that
%this is actually the mirror image of what I should want. Will repair that
%below. 
trodesOmnetics = [...
    8:1:23;
    7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24;
    40:1:55;
    39:-1:32,63:-1:56];
% trodesOmneticsSingle = fliplr(trodesOmneticsSingle); %actually realized I
% had the right orientation all along. XD
%generate the inverse for the flipped headstage. add 
% trodesOmneticsFlip = rot90(rot90(trodesOmneticsSingle))+32;
% %combine the two headstages
% trodesOmnetics = [trodesOmneticsSingle;trodesOmneticsFlip];

%generate master sheet that will hold everything
master = zeros(length(myShankDesig),6);
master(:,1) = myShankDesig;
master(:,2) = neuroTechShankDesig;
%now I want to find the appropriate neuronexus numbers for the same
%positions
for i = 1:64
    posFinder = find(neuroTechDesig==master(i,2));
    master(i,3) = neuraLynxDesig(posFinder);
end
%now convert this to the omnetics connector to get trodes hardware
%designations
for i = 1:64
    posFinder = find(neuralynxOmnetics == master(i,3));
    master(i,4) = trodesOmnetics(posFinder);
end



