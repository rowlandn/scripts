% Nexload.m
% FIRST FILE CALLED BY Nexmerge.m/Spontmerge.m
%
% SHOULD ADD routine for renaming variables
% SHOULD ADD checking if is an EVENT NEX (accelerometer) from OLS using strfind()
% 

TOPDIR = 'G:\Kevinswork\DBS\';

cd( TOPDIR)

if ~exist('n') % controls how many successive LOADS
    n=1;
end

if (~exist('nexpath'));
    [Nexdir(1).name, nexpath] = uigetfile({'*.*',  'All Files (*.*)'});
end

nid=fopen([nexpath Nexdir(n).name],'r');   % Open selected file READ ONLY
if (nid == -1)
    error(sprintf('Cannot open file:\n %s',[nexpath Nexdir(n).name]));
    return
end
% nexfilename = Nexdir(n).name; % use to repeatedly update CURRENT FILE NAME
% if ~isi_batching
%     basefilename = nexfilename; % save for later
% end
    
% *******************************************************************
% READING/LOADING *.nex file FILE HEADER fields (add more GUI / Error Checking)
% *******************************************************************
% NOTE: after READING the file at X precision -> ML converts to DOUBLE or CHAR matrix
Nexfile.MagicNumber = fread(nid, 1, 'int32');
if Nexfile.MagicNumber ~= 827868494
    error('This file is NOT a NEX file!');
    return
end
Nexfile.Version = fread(nid, 1, 'int32');
Nexfile.Comment(1,:) = fread(nid, 256, 'uchar=>char'); % READ as 256 unsigned chars -> char
Nexfile.Frequency = fread(nid, 1, 'float64');
Nexfile.Beg = fread(nid, 1, 'int32');
Nexfile.End = fread(nid, 1, 'int32');
Nexfile.NumVars = fread(nid, 1, 'int32');
Nexfile.NextFileHeader = fread(nid, 1, 'int32');
Nexfile.Padding = fread(nid, [1, 256], 'uchar=>char');
fclose(nid); % EXIT FILE READING (after reading file header)

% USE NEX AUTHOR'S MATLAB FUNCTION TO GRAB: #vars, varnames, vartypes:
% THIS GRABS ALL VARIABLES!! (not just type 0/1)
[nex_nvar nex_names nex_types] = nex_info([nexpath Nexdir(n).name]) % grab NUM vars, NAMES vars, TYPES vars

% *******************************************************************
% LOADING *.nex file VARIABLE HEADER field variables
% *******************************************************************

% NOW BEGIN VECTORIZING (or CELLS) NEX VARIABLE HEADERS:
% Keeps "reopening" the file because of some strange loopspace error...
for i = 1:nex_nvar
    readoffset = 544 + (i-1)*208;
    nid=fopen([nexpath Nexdir(n).name],'r');
    fseek(nid, readoffset, 'bof');
% READ EXACT COPY OF EACH VARIABLE HEADER
    Nexvar.Type(i,:) = fread(nid, 1, 'int32');
    Nexvar.Version(i,:) = fread(nid, 1, 'int32'); % currently seems to be a constant "100"
    Nexvar.Name(i,:) = fread(nid, 64, 'uchar=>char');
    % RLM: CHECK FOR UNIQUE NAME !!!
%%%% Nexvar.DataOffset CALCULATED LATER (INCREMENTED VALUE JUST FOR HEXEDITOR DEBUGGING)
    Nexvar.DataOffset(i,:) = fread(nid, 1, 'int32'); %  [REALLY: index of 1st associated data point in file]
    Nexvar.Count(i,:) = fread(nid, 1, 'int32'); % CHECK THIS FOR ACCURACY LATER
    % BELOW ALL IS FINE FOR NOW:
    % next 6 are usually useless, unless reading in *.nex file with 3-D neuron data
    Nexvar.WireNumber(i,:) = fread(nid, 1, 'int32');
    Nexvar.UnitNumber(i,:) = fread(nid, 1, 'int32');
    Nexvar.Gain(i,:) = fread(nid, 1, 'int32');
    Nexvar.Filter(i,:) = fread(nid, 1, 'int32');
    Nexvar.XPos(i,:) = fread(nid, 1, 'float64');
    Nexvar.YPos(i,:) = fread(nid, 1, 'float64');
    % next 5 are only useful in CONTINUOUS/WAVEFORM/MARKER variables, NOT in neuronal/event
    Nexvar.WFrequency(i,:) = fread(nid, 1, 'float64');
    Nexvar.ADtoMV(i,:) = fread(nid, 1, 'float64');
    Nexvar.NPointsWave(i,:) = fread(nid, 1, 'int32');
    % the next 76 bytes are ALWAYS BLANK (make '0') UNLESS using MARKERS
    Nexvar.NMarkers(i,:) = fread(nid, 1, 'int32');
    Nexvar.MarkerLength(i,:) = fread(nid, 1, 'int32');
    Nexvar.Padding(i,:) = fread(nid, [1, 68], 'uchar=>char');
    fclose(nid); % EXIT FILE READING (after reading all data points)
    % *******************************************************************
    % LOADING the ACTUAL DATA values for Type 0/1 TIMESTAMPS
    % *******************************************************************
    % USE NEX AUTHOR'S MATLAB FUNCTIONS AS HELPERS:
    
    switch Nexvar.Type(i,1)  % Nx1 array of types in file.
    % Nexdata.TS_count = NUMBER OF TIMESTAMPS, Nexdata.TS = TIMESTAMPS in DECIMAL SECONDS (x.xxxxx)
        case {0, 1}   % READ NEURON OR EVENT: one timestamp array
            % NOTE: values of nex_ts have been DIVIDED by Nexfile.Frequency to return SECONDS
            %%% Neuron/Continuous are DIFFERENT sampling, but stored same time base.
            [Nexdata.TS_count(i,1), Nexdata.TS{i,1}] = nex_ts([nexpath Nexdir(n).name], Nexvar.Name(i,:));   % should create Nx64 char arrays
            Nexdata.Freq(i,1) = Nexfile.Frequency;
            if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
            end
            if Nexdata.TS_count(i,1) ~= Nexvar.Count(i,1)
                disp(sprintf('HEADER %d timestamps not equal to ACQUIRED %d timestamps!\n', Nexvar.Count(i,1), Nexdata.TS_count(i,1)));
            end
        case 2  % READ INTERVAL: beginning & end TIMESTAMP ARRAYS
            [Nexdata.TS_count(i,1), Nexdata.TS{i,1}, Nexdata.TS2{i,1}] = nex_int([nexpath Nexdir(n).name], Nexvar.Name(i,:));
            Nexdata.Freq(i,1) = Nexfile.Frequency;
            if Nexdata.TS2{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
            end
        case 3  % READ WAVEFORM: (ts, ad)
            [Nexdata.Freq(i,1), Nexdata.TS_count(i,1), Nexdata.TS{i,1}, Nexdata.FragSizes{i,1}, Nexdata.Cont{i,1}] ...
                = nex_wf([nexpath Nexdir(n).name], Nexvar.Name(i,:));
            % Nexdata.Cont{i,1} IS NF AD/timestamps x N timestamps/waveforms
            if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
            end
        case 4  % READ POPULATION VECTOR
            % SKIP THIS FOR NOW
        case 5  % READ CONTINUOUS A/D: ts, 
            [Nexdata.Freq(i,1), Nexdata.TS_count(i,1), Nexdata.TS{i,1}, Nexdata.FragSizes{i,1}, Nexdata.Cont{i,1}] = nex_cont([nexpath Nexdir(n).name], Nexvar.Name(i,:));
            % NOW CREATE Nexdata.FragIndices from Nexdata.FragSize
            Nexdata.FragIndices{i,1}(1,1) = 0;
            for j = 1:length(Nexdata.FragSizes{i,1})
                Nexdata.FragIndices{i,1}(1,j+1) = Nexdata.FragIndices{i,1}(1,j) + Nexdata.FragSizes{i,1}(1,j);
            end
            if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
            end
        case 6  % READ MARKER
            [Nexdata.TS_count(i,1), Nexdata.NumFields(i,1), Nexdata.FieldLength(i,1), Nexdata.TS{i,1}, Nexdata.FieldNames{i,1}, Nexdata.Cont{i,1}] ...
                = nex_marker([nexpath Nexdir(n).name], Nexvar.Name(i,:));
            % FOR MARKERS: Nexdata.Cont{i,1} contains N x NL x NM array
            % n = num markers x chars/field x fields/marker
            % nm = fields/marker, nl = chars/field
            if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
            end
    end
end
disp('NEX FILE SUCCESSFULLY LOADED')
% ################################################################################
% ################################################################################
% ################################################################################

% RLM: maybe should add a sub-function to deblank/pad varnames,