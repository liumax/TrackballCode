%
%[stimulus] = read_entire_spr_file(sprfile)
%
%
%       FILE NAME       : RT WSTRF DB
%       DESCRIPTION     : Real Time spectro-temporal receptive field
%			  Uses Lee/Schetzen Aproach via Specto-Temporal Envelope
%			  For dB Amplitude Sound distributions 
%
% SpecFile : Spectral Profile File
% sprtype : SPR File Type : 'float' or 'int16'
%         Default=='float'. This function assumes it is a float file.	
%
% stimulus : mxn matrix, which is the entire dmr envelope file.
%
% To get things to match up with Tatyana's code you use the following line
% to compute the STA:
%
% [taxis,faxis,STRF1,PP,Wo1,No1,SPLN] = get_sta_filter(file, 0.005, 0.095, spet, trigger, fs, SPL, MdB, tbins, fbins);
%
% tbins = 20 and fbins = 25 in every case
%
% caa 12/15/06
function [stimulus] = read_entire_spr_file(sprfile, NF, NT)

Sound = 'MR';
ModType = 'dB';
sprtype='float';


%Loading Parameter Data
index = findstr(sprfile,'.spr');
% paramFile = [sprfile(1:index(1)-1) '_param.mat'];
% f = ['load ' paramFile];
% eval(f);
% clear App  MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 f2 
% clear Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn

if nargin == 1
    NF = 64;
    NT = 67;
end
 
%Opening Spectral Profile File
fid = fopen(sprfile);
frewind(fid);
MdB = 40;
RMSP = -MdB/2;
stimulus = [];
while ( ~feof(fid) )
   [s1, count] = fread(fid,NT*NF,'float');
   if count == NF * NT
      s1 = reshape(s1,NF,NT);
      s1 = MdB * s1 - RMSP;
      stimulus = [stimulus s1];
   else
       if mod(count, NF) ~= 0
           remove = mod(count, NF);
           s1(end-remove+1:end) = [];
       end
       s1 = reshape(s1, NF, length(s1)/NF);
       s1 = MdB * s1 - RMSP;
       stimulus = [stimulus s1];

   end % (if)
end % (while)

[nf, ntrials] = size(stimulus);

fprintf('\n#Freqs = %.0f, #trials = %.0f\n\n', nf, ntrials);


%Closing all opened files
fclose('all');
