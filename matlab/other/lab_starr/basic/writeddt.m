function writeddt(neuronal_data, actualsamplerate, filepathandname, varargin)
% function roryddt(neuronal_data, actualsamplerate, fileandpathname)
% CALLED BY: whatever wants to use it! (iopdata.m)
% Rory Michaelis
% OUT: Write unit data in DDT File Format to the "filepathandname" location using input parameters
% IN:   1. neuronal_data = 2 dimensional array, smaller dimension is CHANNELS (typ 1-8), longer is DATA.
%                           CURRENTLY, must provide SCALED DATA! ( rng: +/-32,767, +/-2048 )
%       2. actualsamplerate = the REAL SAMPLERATE used by DAQ or soundcard (may vary)
%       3. filepathandname = datapath,fileroot,&incrementer WITHOUT extension
%       4. varargin{1} = bitspersample (16 if ommitted)
%       5. varargin{2} = filecomment (blank or some default string if ommitted)
% NOTE: could add other inputs for file, channel gain also!
%global data_filepath filename incr % OLD STUFF FROM INTRAOP

% DDT HEADER SECTION: 432 bytes, 
h_version = 102; % int32: 100-102 are possible as of Oct 2003.
h_dataoffset = 432; % int32: assume unit data begins IMMEDIATELY after header bytes
% NOTE: example file "test.ddt" uses 432 byte dataoffset!
h_frequency = actualsamplerate; % double (8 byte): from iopdaq, soundcard or DAQ's actual samplerate
if length(size(neuronal_data)) > 2
    error('writeddt() only writes DDT files for 2 dimensional data arrays.  Not sure how to handle 2+ dimensions!');
    return
end
[numrows numcols] = size(neuronal_data); % int32: currently will only be ONE channel recorded!
if numrows > numcols % successive data is down the rows...
    NPOINTS = numrows;
    h_numchannels = numcols;
    neuronal_data = shiftdim(neuronal_data,1); % reverse row/column for WRITING DATA INTERLEAVED (down col)
else % successive data is across the cols...
    NPOINTS = numcols;
    h_numchannels = numrows;
end
h_clock = fix(clock); % fetch current time in 1x6 array
h_gain = 1; % pre-amp gain for all channels, default is 1
if nargin >= 5
    h_comment = varargin{2};
else
    h_comment = 'Created by abfconv2.m';
end
if nargin >= 4
    h_bitspersample = varargin{1};
else
    h_bitspersample = 16; % 12 or 16 bit resolution (originally DDT just used 12 bits!)
end
h_channelgain(1:h_numchannels) = 1; % 1,2,5,10,20,50, or 100 are legit

% WRITE FILE
fid = fopen([filepathandname '.ddt'], 'w');
fwrite(fid, h_version, 'int32'); % 4 bytes
fwrite(fid, h_dataoffset, 'int32'); % 4 bytes
fwrite(fid, h_frequency, 'double'); % 8 bytes
fwrite(fid, h_numchannels, 'int32'); % 4 bytes
fwrite(fid, h_clock(1), 'int32'); % 4 bytes, yr
fwrite(fid, h_clock(2), 'int32'); % 4 bytes, mon
fwrite(fid, h_clock(3), 'int32'); % 4 bytes, day
fwrite(fid, h_clock(4), 'int32'); % 4 bytes, hr
fwrite(fid, h_clock(5), 'int32'); % 4 bytes, min
fwrite(fid, h_clock(6), 'int32'); % 4 bytes, sec
fwrite(fid, h_gain, 'int32'); % 4 bytes
h_comment(length(h_comment)+1:128) = ' '; % pad 128 byte comment
fwrite(fid, h_comment, 'char'); % 1 byte * 128
fwrite(fid, h_bitspersample, 'uchar'); % 1 byte
h_channelgain(length(h_channelgain)+1:64) = 1; % pad 64 byte channel gain
fwrite(fid, h_channelgain, 'uchar'); % 1 byte * 64
fwrite(fid, zeros(1,191), 'uchar'); % 1 byte * 191 = 191 byte padding
% DDT DATA SECTION
% TRANSPOSE DATA MATRIX SO CHANNELS ARE INTERLEAVED !!! (fwrite() writes DOWN column vectors)
fwrite(fid, neuronal_data, 'int16'); % writes the NPOINTS 16 bit A/D values (interwoven)
fclose(fid);
disp(['FOUND ' num2str(h_numchannels) ' NEURON/ACCEL CHANNEL(s) with ' num2str(NPOINTS) ' DATA POINTS each']);
disp(['Successful DDT file creation at ' filepathandname '.ddt']);
