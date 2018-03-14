%function []=float2sw(infile, outfile, initbuffer, outbuffer)
%
%       FILE NAME       : FLOAT 2 INT 16
%       DESCRIPTION     : Converts a binary 'float' file to binary 'int16'
%                         signed word file
%
%       infile          : Input file name
%       outfile         : Output file name
%       initbuffer      : Number of initial zeros to pre-append to file.
%                         This is used for DVD-Audio .wav files.
%       outbuffer       : Number of zeros to append to the end of the file.
%
function [] = float2sw(infile, outfile, initbuffer, outbuffer)

if ( nargin == 2 )
   initbuffer = 0;
   outbuffer = 0;
end

if ( nargin == 3 )
   outbuffer = 0;
end


%Opening Files
fidin = fopen(infile,'r');
fidout = fopen(outfile,'w');


% Finding Max over the entire file
% for later Normalizing

maxval = eps; % default value

while ~feof(fidin)
   temp = fread(fidin,1024*128,'float');
   maxval = max([maxval abs(temp')]);
end


%Converting to 'int16'

fseek(fidin,0,-1);  % go to beginning of input file

% Write out initial zeros
if ( initbuffer ~= 0 )
   fwrite(fidout, zeros(initbuffer,1), 'int16');
end   

% write out the signed int16 file 
while ~feof(fidin)
   temp = fread(fidin, 1024*128, 'float');
   temp = round(temp ./ maxval * 1024 * 32 * 0.99);
   fwrite(fidout, temp, 'int16');
end

% Write out initial zeros
if ( outbuffer ~= 0 )
   fwrite(fidout, zeros(outbuffer,1), 'int16');
end   

%Closing Files
fclose('all');





