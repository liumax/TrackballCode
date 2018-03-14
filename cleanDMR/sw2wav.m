function sw2wav(chan1file, Fs, nbits, chan2file)
%sw2wav Convert signed integer files to .wav sound files.
%
%   sw2wav(chan1file, Fs, nbits, chan2file) writes data in chan1file.sw
%   and chan2file.sw to a Windows WAVE file specified by the 
%   file name chan1file.wav, with a sample rate of Fs Hz and
%   with nbits precision.
%
%   chan1file, chan2file : a string specifying the file name 
%   of the signed integer files.
%
%   Fs : sampling rate in Hz. The default is 96000 Hz.
%
%   nbits : precision of input files. Default is 16. Other options
%           are 20 and 24.
%
%   sw2wav(chan1file) assumes Fs = 96000 Hz and nbits = 16 and creates 
%   a mono .wav file at 16 bit resolution.
%
%   sw2wav(chan1file, [], [], chan2file) assumes Fs = 96000 Hz
%   and nbits = 16 and creates a stereo .wav file.
%
%   The number of samples in chan1file must equal the number
%   of samples in chan2file or an error is returned.
%
%   For usual usage chan1file is a signal file and chan2file 
%   is a file of triggers corresponding to the signal file.
%
%   See also WAVREAD, AUWRITE.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.14 $  $Date: 2002/04/06 22:40:51 $

%   D. Orofino, 11/95
%   C. Atencio 4/24/03



% Parse inputs:

if ( nargin < 1 & nargin > 4)
   error('Need 1 or 4 input args.');
end

if ( nargin == 1 )
   if ( isstr(chan1file) )
      Fs = 96000;
      channels = 1;
      nbits = 16;
   else
      error('chan1file must be a string.');
   end
elseif ( nargin == 2 )
   if ( isstr(chan1file) )
      channels = 1;
      nbits = 16;
      if ( isempty(Fs) )
         Fs = 96000;
      end
   else
      error('chan1file must be a string.');
   end
elseif ( nargin == 3 )
   if ( isstr(chan1file) )
      channels = 1;
      if ( isempty(Fs) )
         Fs = 96000;
      end
      if ( isempty(nbits) )
         nbits = 16;
      end
   else
      error('chan1file must be a string.');
   end
elseif ( nargin == 4 )
   if ( isstr(chan1file) & isstr(chan2file) )
      channels = 2;
      if ( isempty(Fs) )
         Fs = 96000;
      end
      if ( isempty(nbits) )
         nbits = 16;
      end
   else
      error('chan1file and chan2file must be strings.');
   end
end


% % Specify the precision of the input/output files.
% nbits = 16;

if ( nbits == 16 )
   dtype = 'int16';
elseif ( nbits == 20 )
   dtype = 'bit20';
elseif ( nbits == 24 )
   dtype = 'bit24';
else
   error('nbits must be 16, 20, or 24 only.');
end


% % -------------------- process signal file data and filename -------------
% if ( findstr(chan1file,'.sw') )
%    chan1swfile = chan1file;
%    [chan1fidin, err] = fopen(chan1swfile,'r');
%    error(err);
% else
%    chan1swfile = [chan1file '.sw'];
%    [chan1fidin, err] = fopen(chan1swfile,'r');
%    error(err);
% end



[chan1fidin, err] = fopen(chan1file,'r');
error(err);

% find number of samples in chan1file
samples = 0;
while ~feof(chan1fidin)
   chan1data = fread(chan1fidin, 1024*512, dtype);
   samples = samples + length(chan1data);
end

clear('chan1data');
fclose(chan1fidin);



% --------------------- process stereo data ------------------
if ( nargin == 4 )

%    if ( findstr(chan2file,'.sw') )
%       chan2swfile = chan2file;
%       [chan2fidin, err] = fopen(chan2swfile,'r');
%       error(err);
%    else
%       chan2swfile = [chan2file '.sw'];
%       [chan2fidin, err] = fopen(chan2swfile,'r');
%       error(err);
%    end



   [chan2fidin, err] = fopen(chan2file,'r');
   error(err);

   % find number of samples in chan2file
   chan2samples = 0;
   while ~feof(chan2fidin)
      chan2data = fread(chan2fidin, 1024*512, dtype);
      chan2samples = chan2samples + length(chan2data);
   end

   clear('chan2data');
   fclose(chan2fidin);

   if ( samples ~= chan2samples )
      error('chan1file and chan2file must have same number of samples.');
   end

end



% [samples, channels] = size(y);


% Determine number of bytes in chunks
% (not including pad bytes, if needed):
% ----------------------------------
%  'RIFF'           4 bytes
%  size             4 bytes (ulong)
%  'WAVE'           4 bytes
%  'fmt '           4 bytes
%  size             4 bytes (ulong)
% <wave-format>     14 bytes
% <format_specific> 2 bytes (PCM)
%  'data'           4 bytes
%  size             4 bytes (ulong)
% <wave-data>       N bytes
% ----------------------------------

bytes_per_sample = ceil(nbits/8);
total_samples    = samples * channels;
total_bytes      = total_samples * bytes_per_sample;

riff_cksize = 36+total_bytes;   % Don't include 'RIFF' or its size field
fmt_cksize  = 16;               % Don't include 'fmt ' or its size field
data_cksize = total_bytes;      % Don't include 'data' or its size field

% Determine pad bytes:
data_pad    = rem(data_cksize,2);
riff_cksize = riff_cksize + data_pad; % + fmt_pad, always 0

% Open file for output:
[fidout, err] = OpenWaveWrite(chan1file);
error(err);

% Prepare basic chunk structure fields:
ck = []; 
ck.fid = fidout; 
ck.filename = chan1file;

% Write RIFF chunk:
ck.ID   = 'RIFF';
ck.Size = riff_cksize;
error(write_ckinfo(ck));

% Write WAVE subchunk:
ck.ID   = 'WAVE';
ck.Size = [];  % Indicate a subchunk (no chunk size)
error(write_ckinfo(ck));

% Write <fmt-ck>:
ck.ID   = 'fmt ';
ck.Size = fmt_cksize;
error(write_ckinfo(ck));

% Write <wave-format>:
fmt.filename        = chan1file; % this file name will be used to create the .wav output file

if nbits == 32,
    fmt.wFormatTag  = 3;            % Data encoding format (1=PCM, 3=Type 3 32-bit)
else
    fmt.wFormatTag  = 1;            
end

fmt.nChannels       = channels;     % Number of channels
fmt.nSamplesPerSec  = Fs;           % Samples per second
fmt.nAvgBytesPerSec = channels*bytes_per_sample*Fs; % Avg transfer rate
fmt.nBlockAlign     = channels*bytes_per_sample;    % Block alignment
fmt.nBitsPerSample  = nbits;        % standard <PCM-format-specific> info
error(write_wavefmt(fidout,fmt));

% Write <data-ck>:
ck.ID   = 'data';
ck.Size = data_cksize;
error(write_ckinfo(ck));



% ----------------- write out the .wav file -------------------

if ( nargin < 4 ) % process monaural data

   [chan1fidin, err] = fopen(chan1file,'r');
   error(err);

   % write the new .wav file
   samples = 0;
   while ~feof(chan1fidin)
      y = fread(chan1fidin, 1024*256, dtype);
      error(write_wavedat(fidout,fmt,y));   % write <wave-data>, and its pad byte if needed
   end

elseif ( nargin == 4 ) % Process stereo data: two input files

   [chan1fidin, err] = fopen(chan1file,'r');
   error(err);
   [chan2fidin, err] = fopen(chan2file,'r');
   error(err);

   % write the new .wav file
   while ~feof(chan1fidin)
      chan1data = fread(chan1fidin, 1024*256, dtype);
      chan2data = fread(chan2fidin, 1024*256, dtype);
      y = [chan1data(:) chan2data(:)];
      error(write_wavedat(fidout,fmt,y));   % write <wave-data>, and its pad byte if needed
   end

end

   
% Determine if a pad-byte must be appended to data chunk:
if rem(total_bytes, 2) ~= 0,
   fwrite(fidout,0,'uchar');  % write out the one byte
end


% Close files:
fclose('all');

% end of sw2wav()


% ------------------------------------------------------------------------
% Private functions:
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
function [fid,err] = OpenWaveWrite(wavefile)
% OpenWaveWrite
%   Open WAV file for writing.
%   If filename does not contain an extension, add ".wav"

fid = [];
err = '';

if ~isstr(wavefile),
   err = 'Wave file name must be a string.'; return;
end

inddot = findstr(wavefile,'.');

if isempty( inddot )
   wavefile = [wavefile '-stereo' '.wav'];
else
   prefix = wavefile(1:inddot-1);
   wavefile = [prefix '-stereo' '.wav'];
end

% Open file, little-endian:
[fid,err] = fopen(wavefile,'wb','l');

return


% ------------------------------------------------------------------------
function err = write_ckinfo(ck)
% WRITE_CKINFO: Writes next RIFF chunk, but not the chunk data.
%   Assumes the following fields in ck:
%         .fid   File ID to an open file
%         .ID    4-character string chunk identifier
%         .Size  Size of chunk (empty if subchunk)
%
%
%   Expects an open FID pointing to first byte of chunk header,
%   and a chunk structure.
%   ck.fid, ck.ID, ck.Size, ck.Data

errmsg = ['Failed to write ' ck.ID ' chunk to WAVE file: ' ck.filename];
err    = '';

if (fwrite(ck.fid, ck.ID, 'char') ~= 4),
   err = errmsg; return;
end

if ~isempty(ck.Size),
  % Write chunk size:
  if (fwrite(ck.fid, ck.Size, 'ulong') ~= 1),
     err = errmsg; return;
  end
end

return


% ------------------------------------------------------------------------
function err = write_wavefmt(fid, fmt)
% WRITE_WAVEFMT: Write WAVE format chunk.
%   Assumes fid points to the wave-format subchunk.
%   Requires chunk structure to be passed, indicating
%   the length of the chunk.

errmsg = ['Failed to write WAVE format chunk to file' fmt.filename];
err    = '';

% Create <wave-format> data:
if (fwrite(fid, fmt.wFormatTag,      'ushort') ~= 1) | ...
   (fwrite(fid, fmt.nChannels,       'ushort') ~= 1) | ...
   (fwrite(fid, fmt.nSamplesPerSec,  'ulong' ) ~= 1) | ...
   (fwrite(fid, fmt.nAvgBytesPerSec, 'ulong' ) ~= 1) | ...
   (fwrite(fid, fmt.nBlockAlign,     'ushort') ~= 1),
   err=errmsg; return;
end

% Write format-specific info:
if fmt.wFormatTag == 1 | fmt.wFormatTag == 3,
  % Write standard <PCM-format-specific> info:
  if (fwrite(fid, fmt.nBitsPerSample, 'ushort') ~= 1),
     err = errmsg; return;
  end
  
else
  err = 'Unknown data format.';
end

return


% -----------------------------------------------------------------------
function y = PCM_Quantize(x, fmt)
% PCM_Quantize:
%   Scale and quantize input data, from [-1, +1] range to
%   either an 8-, 16-, or 24-bit data range.

% Clip data to normalized range [-1,+1]:
ClipMsg  = ['Data clipped during write to file:' fmt.filename];
ClipWarn = 0;

% Determine slope (m) and bias (b) for data scaling:
nbits = fmt.nBitsPerSample;
m = 2.^(nbits-1);

switch nbits
case 8,
   b = 128;
case {16,24},
   b = 0;
otherwise,
   error('Invalid number of bits specified.');
end

y = round(m .* x + b);

% Determine quantized data limits, based on the
% presumed input data limits of [-1, +1]:
ylim = [-1 +1];
qlim = m * ylim + b;
qlim(2) = qlim(2)-1;

% Clip data to quantizer limits:
i = find(y < qlim(1));
if ~isempty(i),
   warning(ClipMsg); ClipWarn=1;
   y(i) = qlim(1);
end

i = find(y > qlim(2));
if ~isempty(i),
   if ~ClipWarn, warning(ClipMsg); end
   y(i) = qlim(2);
end

return


% -----------------------------------------------------------------------
function err = write_wavedat(fid,fmt,data)
% WRITE_WAVEDAT: Write WAVE data chunk
%   Assumes fid points to the wave-data chunk
%   Requires <wave-format> structure to be passed.

err = '';

if fmt.wFormatTag==1 | fmt.wFormatTag==3,
   % PCM Format
   
   % 32-bit Type 3 is normalized, so no scaling needed.
%    if fmt.nBitsPerSample ~= 32,
%        data = PCM_Quantize(data, fmt);
%    end
   
   switch fmt.nBitsPerSample
   case 8,
      dtype = 'uchar'; % unsigned 8-bit
   case 16,
      dtype = 'short'; % signed 16-bit
   case 24,
	  dtype = 'bit24'; % signed 24-bit
   case 32,
      dtype = 'float'; % normalized 32-bit floating point
   otherwise,
      err = 'Invalid number of bits specified.'; 
      return;
   end
   
   % Write data, one row at a time (one sample from each channel):
   [samples,channels] = size(data);
   total_samples = samples*channels;
   
   if (fwrite(fid, reshape(data',total_samples,1), dtype) ~= total_samples),
      err = 'Failed to write PCM data samples.'; 
      return;
   end
   
%    % Determine # bytes/sample - format requires rounding
%    %  to next integer number of bytes:
%    BytesPerSample = ceil(fmt.nBitsPerSample/8);
%    
%    % Determine if a pad-byte must be appended to data chunk:
%    if rem(total_samples*BytesPerSample, 2) ~= 0,
%       fwrite(fid,0,'uchar');
%    end
   
else
  % Unknown wave-format for data.
  err = 'Unsupported data format.';
end

return

% end of wavwrite.m
