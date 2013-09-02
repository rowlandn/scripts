
% nexinit.m
% -------------------------------------------------------------------
% INITIALIZATION OF NEX FILE STRUCTURE: 
% Header Components (Matlab variables) are listed in physical file order
% Header Components are assumed integral UNLESS stated otherwise in comment
% -------------------------------------------------------------------
% INITIALIZE FILE HEADER STRUCTURE COMPONENTS (544 bytes) = "0220" bytes
% -------------------------------------------------------------------

Nexfile.MagicNumber = ''; % 4 bytes... SPECIAL IDENTIFIER 'NEX1' or int32: 827,868,494
Nexfile.Version = [];     % 4 bytes... 101 in Rob's files. (100 to 104 apparently valid)
Nexfile.Comment = '';     % CHAR 256 bytes... Notes
Nexfile.Frequency = [];   % FPOINT-64b 8 bytes... 20,000hz = hex 40D388000000000000
                    % Seems to be used for ALL TYPES
Nexfile.Beg = [];         % 4 bytes... USUALLY 0
Nexfile.End = [];         % 4 bytes... ~MAX Timestamp + 1 (or could do ffffffff)
Nexfile.NumVars = [];     % 4 bytes... Number of variables (hence variable headers) in file 
% the next 260 bytes are ~ALWAYS BLANK (make '0')
Nexfile.NextFileHeader = [];   % 4 bytes... currently an UNUSED FIELD
Nexfile.Padding = '';     % CHAR 256 bytes... Expansion space FOR FUTURE USE

% -------------------------------------------------------------------
% INITIALIZE VARIABLE HEADER STRUCTURE COMPONENTS (208 bytes each) = "D0" bytes
% Variables are vectorized, each ROW contains a value for a different header
% -------------------------------------------------------------------

Nexvar.Type = [];        % 4 bytes...   0=NEURON, 1=event, 2=interval, 3=waveform
                      % 4=population vector, 5=CONTINUOUS VARIABLE, 6=MARKER
Nexvar.Version = [];     % 4 bytes...  100 in Rob's files. (KEEP CONSTANT I THINK)
% Nexvar.Name may be a PROBLEM.. Maybe undefine it and switch to variable length CELL ARRAY ???
Nexvar.Name = '';        % CHAR 64 bytes... Channel descriptions
Nexvar.DataOffset = [];  % 4 bytes... BYTE # WHERE the Data array for variable BEGINS in file
Nexvar.Count = [];       % 4 bytes... NUMBER of TS (hence events, intervals, fragments (for wf / cont) or weights)
% the next 32 bytes are always blanked out (make '0') UNLESS RARE 3-D NEURON DATA
Nexvar.WireNumber = [];  % 4 bytes... NEURON ONLY, currently UNUSED FIELD
Nexvar.UnitNumber = [];  % 4 bytes... NEURON ONLY, currently UNUSED FIELD
Nexvar.Gain = [];        % 4 bytes... NEURON ONLY, currently UNUSED FIELD
Nexvar.Filter = [];      % 4 bytes... NEURON ONLY, currently UNUSED FIELD
Nexvar.XPos = [];        % FPOINT-64b 8 bytes... NEURON ONLY, electrode (0.0-100.0) used in 3-D
Nexvar.YPos = [];        % FPOINT-64b 8 bytes... NEURON ONLY, electrode (0.0-100.0) used in 3-D
% the next 20 bytes (3 components) are set ONLY IF using CONTINUOUS or WAVEFORM variable
Nexvar.WFrequency = [];  % FPOINT-64b 8 bytes... BOTH, cont sample freq: 1000 = 408F400000000000
Nexvar.ADtoMV = [];      % FPOINT-64b 8 bytes... BOTH, conversion factor: mV=1 or 0.001=in volts
Nexvar.NPointsWave = []; % 4 bytes... BOTH, number of points per wave OR TOTAL POINTS in Continuous
% the next 76 bytes are ALWAYS BLANK (make '0') UNLESS using MARKERS
Nexvar.NMarkers = [];    % 4 bytes... number of FIELDS associated with each MARKER
Nexvar.MarkerLength = [];% 4 bytes... number of characters in each MARKER VALUE
Nexvar.Padding = '';     % CHAR 68 bytes... Expansion space FOR FUTURE USE

% USES following parameters, create if necessary!:
current_path = cd;  % just save current path for the hell of it
if ~exist('merge_nexfiles')
    merge_nexfiles = 0; % DEFAULT, if no instructions given.
end
if ~exist('isi_batching')
    isi_batching = 0; % DEFAULT, if no instructions given.
end

% *******************************************************************
% OPEN *.nex file containing neuron data & probably other data
% *******************************************************************
if isi_batching   % batch load all *.nex files in directory
    loadstring = 'Choose the DIRECTORY containing the *.NEX files used for batch ISI statistics.';
    loadstring = [loadstring '  (Output data file''s name will be based upon directory name)'];
    nexpath = uigetdir(cd, loadstring); % NOTE: no fileseperator '\' is appended 
    cd(nexpath); % change current directory to one containing nex files, modified by Sho (8/1/2007)
    if nexpath == 0
        error('You must select a valid directory!');
    end
    fnd = findstr(filesep, nexpath);
    basefilename = [nexpath(fnd(end)+1:end) '.nex']; % filename = parent DIRNAME.nex
    nexpath = [nexpath filesep]; % NOW APPEND FILE SEPARATOR!
    Nexdir = dir([nexpath '*.nex']);
    NUMFILES = length(Nexdir);
    if NUMFILES == 0
        warning('No *.NEX files found in: %s',nexpath);
        return
    end
%elseif merge_batching  % for non-isi BATCH case.. joins continuous data a zillion times!
% 
else
    loadstring = 'Select the NEX file you would like to work with';
    [nexfilename, nexpath] = uigetfile('*.nex', loadstring);
    if nexfilename == 0 & nexpath == 0
        warning(sprintf('Sorry, you must select a valid NEX file!\n EXITING...'));
        merge_nexfiles = 0;
        is_event = 0;
        break;
    end
    clear Nexdir % make sure no residual structure
    Nexdir(1).name = nexfilename;
    NUMFILES = 1;
end
