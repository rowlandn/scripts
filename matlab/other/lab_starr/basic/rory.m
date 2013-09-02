% Rory L. Michaelis
% 2002.08.26 (2002.08.19)
%
% SCRIPT for joining NEX Neuronal Spike Data Files (*.nex) ...
% ...with Continuously Sampled Data loaded into the Matlab workspace
% The Cont Sampled Data comes from the "abfconv.m" script by Rob Turner
% 
% CURRENTLY MONOLITHIC STYLE (modularity later if necessary)
% NOTE: if examing file via HEX EDITOR: x86 uses LITTLE ENDIAN (01 00 = 1, 00 01 = 256)
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% -------------------------------------------------------------------
% INITIALIZATION OF NEX FILE STRUCTURE: 
% all components (aka variables) are listed in physical order
% all variables are assumed integral UNLESS stated otherwise in comments
% -------------------------------------------------------------------
% INITIALIZE FILE HEADER STRUCTURE COMPONENTS (544 bytes) = 0220 bytes
% -------------------------------------------------------------------

F_MagicNumber = ''; % 4 bytes... SPECIAL IDENTIFIER 'NEX1' is int: 827,868,494
F_Version = '';     % 4 bytes... 101 in Rob's files. (100 to 104 valid)
F_Comment = '';     % CHAR will be 256 bytes
F_Frequency = '';   % FPOINT 8 bytes... 20,000hz = hex 40D388000000000000
F_Beg = '';         % 4 bytes... USUALLY 0
F_End = '';         % 4 bytes... ~MAX Timestamp + 1 (or could do ffffffff)
F_NumVars = '';     % 4 bytes... Number of variables (& variable headers) in NEX file
% below is 260 bytes to ALWAYS blanked out (make '0')
F_NextFileHeader = '';   % 4 bytes... currently UNUSED FIELD
F_Padding = '';     % CHAR will be 256 bytes... FOR FUTURE USE

% -------------------------------------------------------------------
% INITIALIZE VARIABLE HEADER STRUCTURE COMPONENTS (208 bytes) = D0 bytes
% -------------------------------------------------------------------
% These should be vectorized with each ROW containing a different Var Header

V_Type = '';        % 4 bytes...   0=neuron, 1=event, 2=interval, 3=waveform
                      % 4=population vector, 5=continuous variable, 6=marker
V_Version = '';     % 4 bytes...  100 in Rob's files. (KEEP CONSTANT I THINK)
% V_Name is a PROBLEM.. Maybe undefine it and switch to variable length CELL ARRAY ???
V_Name = '';        % CHAR will be  64 bytes... Describes channel
V_DataOffset = '';  % 4 bytes... BYTE # WHERE the Data array for variable BEGINS in file
V_Count = '';       % 4 bytes... NUMBER of TS/events, intervals, waveforms, fragments, or weights
% below is 32 bytes always blanked out (make '0') UNLESS RARE 3-D NEURON DATA
V_WireNumber = '';  % 4 bytes... NEURON ONLY, currently UNUSED FIELD
V_UnitNumber = '';  % 4 bytes... NEURON ONLY, currently UNUSED FIELD
V_Gain = '';        % 4 bytes... NEURON ONLY, currently UNUSED FIELD
V_Filter = '';      % 4 bytes... NEURON ONLY, currently UNUSED FIELD
V_XPos = '';        % FPOINT 8 bytes... NEURON ONLY, electrode (0-100) used in 3-D
V_YPos = '';        % FPOINT 8 bytes... NEURON ONLY, electrode (0-100) used in 3-D
% below is useful stuff ONLY IF using CONTINUOUS or WAVEFORM variable
V_WFrequency = '';  % FPOINT 8 bytes... BOTH, cont sample freq: 1000 = 408F400000000000
V_ADtoMV = '';      % FPOINT 8 bytes... BOTH, Millivolt conversion factor: 1 or 0.001=in volts 
V_NPointsWave = ''; % 4 bytes... WAVEFORM ONLY, number of points per wave
% below is 76 bytes always blanked out (make '0') UNLESS using MARKERS
V_NMarkers = '';    % 4 bytes... number of values associated with each MARKER
V_MarkerLength = '';% 4 bytes... number of characters in each MARKER VALUE
V_Padding = '';     % CHAR will be 68 bytes... FOR FUTURE USE (Should it be zeroed??)

% *******************************************************************
% *******************************************************************
% READING/LOADING *.nex file FILE HEADER fields (add more GUI / Error Checking)
% *******************************************************************

global nexfilename nexpathname;
[nexfilename, nexpathname] = uigetfile('*.nex', 'Select NEX file with ONLY Neuronal Data');
fid=fopen(nexfilename,'r');     % Open selected file READ ONLY
if (fid == -1)
	disp('ERROR:  Cannot open file!');
   return
end
    % NOTE: after reading the file at X precision -> converted to DOUBLE or CHAR matrix  
    F_MagicNumber = fread(fid, 1, 'int32');
    F_Version = fread(fid, 1, 'int32');
    F_Comment(1,:) = fread(fid, 256, '256*uchar=>char');
    F_Frequency = fread(fid, 1, 'float64');
    F_Beg = fread(fid, 1, 'int32');
    F_End = fread(fid, 1, 'int32');
    F_NumVars = fread(fid, 1, 'int32');
    F_NextFileHeader = fread(fid, 1, 'int32');
    F_Padding = zeros(1,256);
fclose(fid);

% *******************************************************************
% LOADING *.nex file VARIABLE HEADER field variables
% (currently works ONLY with *.nex files containing ONLY Type 0/1 NEURON/EVENT DATA CHANNELS)
% *******************************************************************

% USE NEX AUTHOR'S MATLAB FUNCTION TO GRAB: #vars, varnames, vartypes:
%%% EVENTUALLY MAKE THIS INDEPENDENT OF HIS FILES !!!
% THIS ACTUALLY GRABS ALL VARIABLES!! (not just type 0/1)
[nex_nvar nex_names nex_types] = nex_info(nexfilename)

% NOW BEGIN VECTORIZING NEURONAL VARIABLE HEADERS:
%%%% THIS WILL BE BROKEN, IF THE NEX FILE HAS NON-NEURONAL VARIABLES !!!!
V_Type = [V_Type; nex_types];  % n_types   % CONCATENATE nex_types into V_Type header component
V_Name = [V_Name; nex_names];  % n_names   % CONCATENATE nex_names into V_Name header component

for i = 1:nex_nvar  % initialize "nex_var" rows for each neuronal variable header component
    V_Version = [V_Version; 100];      % seems to be constant...
%%%% V_DataOffset is COMPLICATED (INCREMENTED INDEX VALUE IS JUST FOR HEXEDITOR DEBUGGING)
    V_DataOffset = [V_DataOffset; i]  %  [REALLY: index of 1st associated data point in file]
    % ALL IS FINE BELOW:
    % next 6 are typically useless, unless we are reading in *.nex file containing 3-D data already
    V_WireNumber = [V_WireNumber; 0];
    V_UnitNumber = [V_UnitNumber; 0];
    V_Gain = [V_Gain; 0];
    V_Filter = [V_Filter; 0];
    V_XPos = [V_XPos; 0];        % used only in 3-D Neuron
    V_YPos = [V_YPos; 0];        % used only in 3-D Neuron
    % next 5 are only useful in CONTINUOUS/WAVEFORM/MARKER variables, NOT in neuronal
    V_WFrequency = [V_WFrequency; 0];
    V_ADtoMV = [V_ADtoMV; 0];
    V_NPointsWave = [V_NPointsWave; 0];
    V_NMarkers = [V_NMarkers; 0];
    V_MarkerLength = [V_MarkerLength; 0];
% *******************************************************************
% LOADING the ACTUAL data values for Type 0/1 Continuous Data
% *******************************************************************
% USE NEX AUTHOR'S MATLAB FUNCTION TO GRAB TIME STAMPS:
%%% EVENTUALLY MAKE THIS INDEPENDENT OF HIS FILES !!!
    % D_neuron_ts_number = NUMBER OF TIMESTAMPS, D_neuron_ts = TIMESTAMPS in DECIMAL SECONDS (x.xxxxx)
    % Returned values of nex_ts have been DIVIDED by F_Frequency
%%% SHOULD ALL types of TS be unified in one (cell) variable??  Neuron/Continuous are DIFFERENT sampling.
    %% Next line will NOT work if V_Name rows are varying lengths (non-64)...
    [D_neuron_ts_number(i,1) D_neuron_ts(i,:)] = nex_ts(nexfilename, V_Name(i,:));    % should create Nx64 char arrays
    V_Count = [V_Count; D_neuron_ts_number(i,1)];    % V_Count(i) gets the NUMBER of timestamps for each Neuronal Channel
end
V_Padding = zeros(1,68);        % no reason to VECTORIZE padding!
    
% HERE IS WHERE WE SHOULD ADD "NEX FILE" DATA LOGIC FOR TYPE 2,3,4,5,6 variables !!

% *******************************************************************
% LOADING *.mat file variables as Type 3 WAVEFORM DATA CHANNELS 
% *******************************************************************

% This was automated to find the nexfilename's corresponding matfilename.
load(strcat(nexfilename(1:length(nexfilename)-3), 'mat'));
%[matfilename, matpathname] = uigetfile('*.mat', 'Select MAT file with Continuously Sampled Data');
%load(matfilename);          % GRABS CONTINUOUS DATA: n_emg, emg_chan, time     

for cont_index = 1:n_emg  % n_emg (1x1) = number of continuous sampled channels in *.mat file
    V_Type = [V_Type; 3];           % OR should this be 3 (waveform) ??
    V_Version = [V_Version; 100];   % seems to be constant...
%%% how should cont channels be NAMED?? it's not stored in *.mat file !! % PROMPT FOR A NAME ???
    V_Name = char(V_Name(:,:), strcat('emg_chan', num2str(cont_index))); % use "char()" to add trailing spaces
%%% V_DataOffset is COMPLICATED (INCREMENTING VALUE IS JUST FOR HEXEDITOR DEBUGGING)
    V_DataOffset = [V_DataOffset; cont_index+F_NumVars];  % [REALLY: index of 1st associated data point in file]
%%% V_Count should be made independent, based upon each separate READ / DATA CHANNEL INPUT
    V_Count = [V_Count; length(emg_chan(cont_index,:))];  % typically ALL CONTINUOUS VARIABLES SAME LENGTH
    % next 6 are USELESS (only for neuronal data)
    V_WireNumber = [V_WireNumber; 0];
    V_UnitNumber = [V_UnitNumber; 0];
    V_Gain = [V_Gain; 0];
    V_Filter = [V_Filter; 0];
    V_XPos = [V_XPos; 0];
    V_YPos = [V_YPos; 0];
    % next 5 are only useful in CONTINUOUS/WAVEFORM/MARKER variables, NOT in neuronal
    V_WFrequency = [V_WFrequency; 1000];    % Get user input ??????
% THIS IS NOT TO MILLIVOLTS!! I USED x100 scale factor WRITING VALS cuz integral A/D values!
    V_ADtoMV = [V_ADtoMV; 0.005];             % Get user input ??????
    V_NPointsWave = [V_NPointsWave; 1];     % ONLY FOR WAVEFORM? (type 3)
    V_NMarkers = [V_NMarkers; 0];
    V_MarkerLength = [V_MarkerLength; 0];        
% *******************************************************************
% LOADING the ACTUAL data values for Type 3/5 Continuous Data
% *******************************************************************
%%% Currently REQUIRES each Continuous Channel have the SAME number of TS
%%% ASSUMES 1 timestamp / waveform, as that is the format of Rob's Matlab datastructures!
    D_cont_ts(cont_index,:) = time;  % time (1xX) = X time stamps in DECIMAL SECONDS (x.xxxxx):
    D_cont_wf(cont_index,:) = emg_chan(cont_index,:);  % emg_chan (nxX) = N channels with X values each channel:
%    D_neuron_ts(cont_index+F_NumVars,:) = D_cont_ts(cont_index,:)  % THIS WILL BE THE REAL TS ARRAY (for future use)
end
D_neuron_number=F_NumVars;
F_NumVars = F_NumVars + n_emg;   % set CORRECT number of variables for the file...

% *******************************************************************
% CALCULATE/LOAD the OFFSET FIELDS for 0/1 & Type3/5 Continuous Data
V_DataOffset(1,:) = 544 + F_NumVars*208; % initial offset, just after all F_NumVars headers
%%% SHOULD "V_Count(i,:)" be used for offsets/writing ????
for i = 2:D_neuron_number+1   % assume the FIRST D_neuron_number of V_Type are "0" or "1"
%%% PERHAPS USE A CASE SWITCH STATEMENT FOR THE SIX TYPES OF VARIABLES?
% WE MUST ONLY WORRY ABOUT 0/1 and 3 for now.. 5 soon.. and 2 is EASY.. 4/6 ?????
        V_DataOffset(i,:) = V_DataOffset(i-1,:)+4*D_neuron_ts_number
end
%%% WILL WE EVER CALCULATE:  NPointsWave ????
for i = D_neuron_number+2:F_NumVars
    switch V_Type(i,:)  % continue scanning through the V_type vector
        case 3
        V_DataOffset(i,:) = V_DataOffset(i-1,:) + 4*length(D_cont_ts(i-1,:)) + ...
            2*length(D_cont_wf(i-1,:))
        case 5
        disp('incomplete')
    end
end

% ###################################################################
% ###################################################################
% WRITING the FILE HEADER crap to construct a new NEX file (this is perfect)
% ###################################################################
%%%%%%%%%%%%%%%%%%%%%%%%%% Make nice file SAVE dialogue later!
fid = fopen('crap.nex', 'r+');   % for now, the file must be PRE-EXISTING
fwrite(fid, F_MagicNumber, 'int32');
fwrite(fid, F_Version, 'int32');
%%% SPACES SHOULD NOT MATTER, COMMENT IS VERBATIM FROM THE *.nex file
%%% DISPLAY CURRENT CHANNEL NAMES ???  (AT BOTTOM prompt to edit them maybe...)
disp(sprintf('\n\n===================================================================================='))
disp(sprintf('\nNow preparing to write NEW nex file...\n\nThe current CHANNEL NAMES are:\n'))
disp(char(V_Name))
disp(sprintf('\nThe current NEX FILE COMMENT is:\n%s\n\n', deblank(F_Comment)))
choice=input('WOULD YOU LIKE TO CREATE A NEW NEX FILE COMMENT???\n If so, type "yes" and press ENTER. If not, just press ENTER\n # > ', 's')
if strcmp(choice, 'yes')
    F_Comment = input('\n\nType new nex file comment now.  Feel free to copy and paste. \n REMEMBER, ONLY 256 CHARACTERS WILL BE WRITTEN!!!\n # > ', 's')
end
commentsize=length(F_Comment)
if commentsize < 256
    empty = 256 - commentsize       % Should "empty" be a zeroed matrix instead?
    for i = 1:commentsize;   % does this work going from 1->0 if commentsize=0 ?
        fwrite(fid, F_Comment(i),'uchar');
    end
    empty = zeros(1, empty);     % THIS SHOULD BE FASTER THAN A LOOP WRITE!!!
    fwrite(fid, empty, 'uchar');
else 
    fwrite(fid, F_Comment(1:256), 'uchar'); % if Comment GREATER than 256, truncate it
end
% back to simplicity
fwrite(fid, F_Frequency, 'float64');
fwrite(fid, F_Beg, 'int32');
fwrite(fid, F_End, 'int32');         % THIS IS THE ONLY FILE HEADER VALUE TO FIX !!!
fwrite(fid, F_NumVars, 'int32');     % THIS ONE CAN BE UPDATED AS WELL
fwrite(fid, F_NextFileHeader, 'int32');
fwrite(fid, F_Padding, 'uchar');

% ###################################################################
% WRITING the basic VARIABLE HEADER crap to construct a new NEX file
% ###################################################################
% control this by a loop counter variable (matrix row = channel variable)
for V_index = 1:F_NumVars
    fwrite(fid, V_Type(V_index), 'int32');
    fwrite(fid, V_Version(V_index), 'int32');
%%% code to WRITE/ALTER variable length channel names: DEBLANK the padded strings!
disp(sprintf('Variable/Channel Number:  %d  is named:\n%s\n', V_index, deblank(V_Name(V_index,:))))
choice=input('WOULD YOU LIKE TO CHANGE THIS???\n If so, type "yes" and press ENTER. If not, just press ENTER\n # > ', 's')
%%% SHOULD THE NEW NAME BE SAVED BACK INTO THE V_Name() DATA STRUCTURE ???
%%% MUST ALL BE SAME # OF COLUMNS !!!
if strcmp(choice, 'yes')
    varname = input('\nType new VARIABLE NAME now.  Feel free to copy and paste. \n REMEMBER, ONLY 64 CHARACTERS WILL BE WRITTEN!!!\n # > ', 's')
else 
    varname = V_Name(V_index,:);
end
    varname=deblank(varname);
    varnamesize=length(varname)
    if varnamesize < 64
        empty = 64 - varnamesize
        for i = 1:varnamesize          % does this work going from 1->0 if empty name (namesize=0) ?
            fwrite(fid, varname(i), 'uchar');
        end
        empty = zeros(1, empty);     % THIS SHOULD BE FASTER THAN A LOOP WRITE!!!
        fwrite(fid, empty, 'uchar');
    else 
        fwrite(fid, V_Name(V_index,1:64), 'uchar'); % if Comment GREATER than 64, truncate it
    end
    % back to simplicity
    fwrite(fid, V_DataOffset(V_index), 'int32');
    fwrite(fid, V_Count(V_index), 'int32');
    fwrite(fid, V_WireNumber(V_index), 'int32');
    fwrite(fid, V_UnitNumber(V_index), 'int32');
    fwrite(fid, V_Gain(V_index), 'int32');
    fwrite(fid, V_Filter(V_index), 'int32');
    fwrite(fid, V_XPos(V_index), 'float64');
    fwrite(fid, V_YPos(V_index), 'float64');
    fwrite(fid, V_WFrequency(V_index), 'float64');
    fwrite(fid, V_ADtoMV(V_index), 'float64');
    fwrite(fid, V_NPointsWave(V_index), 'int32');
    fwrite(fid, V_NMarkers(V_index), 'int32');
    fwrite(fid, V_MarkerLength(V_index), 'int32');
    fwrite(fid, V_Padding, 'uchar');
end

% ###################################################################
% WRITING the VARIABLE DATA to construct a new NEX file
% ###################################################################

% NOW WRITE THE NEURONAL DATA (SHOULD MULT BY F_Frequency)
% It may be a good idea to just ALTER nex_ts() so it does NOT divide the raw Timestamp Values
for i=1:D_neuron_number
    fwrite(fid, F_Frequency*D_neuron_ts(D_neuron_number,:).','int32');
end
% NOW WRITE THE WAVEFORM/CONTINUOUS DATA
end
for i=D_neuron_number+1:F_NumVars
    fwrite(fid, F_Frequency*D_cont_ts(i-D_neuron_number,:).', 'uint32');
%%% Maybe later WF will have to be INSIDE a cell array, holding the wf array!
%%% LATER ADD DYNAMIC SCALING !!!
    fwrite(fid, 200*D_cont_wf(i-D_neuron_number,:).', 'uint16');
end

fclose(fid);

% ##########################################################################
% *******************************************************************
%           JUST SOME NOTES TO MYSELF:
% *******************************************************************
% EXAMPLES
% *******************************************************************
% DEBUG STATEMENT:  fprintf('comment size %d\n', commentsize)

% concatenate array:   % Type = [Type; VAL]
% decimate array:      % Type(1,:) = []
% CONVERT NUMBER VALUE TO TEXT: num2str(version)
% NOTE:  when WRITING variables to NEX file, use these precisions:
% 1byte -> 'uchar' 4bytes -> 'int32' 8bytes -> 'float64'
% disp(char(MagicNumber))

% INPUT / OUTPUT FUNCTIONS:
% OUTPUT:
% save (simple), fprintf (versatile), dlmwrite, diary (script copies screen output)
% fwrite(fid, variable, precision) %fseek(  %fscan( 
%    STATUS = FSEEK(FID, OFFSET, ORIGIN) repositions the file position
%    ORIGIN values are interpreted   
% BEGIN: 'bof' or -1  CURRENT: 'cof' or  0  END: 'eof' or  1

% NeuroExplorer Executable: IMPORT FORMAT INFO:  
%   NEURON Text COLUMN:  NAMEa\t NAMEb\n TS1a\t TS1b\n TS2a\tTS2b\n
%   OR  LIST: NAMEa\t TS1a \n NAMEb\t TS1b \n NAMEb\t  TS2b \n NAMEb\t TS3b \n
%   CONTINUOUS Text: NAME\n VAL1\n VAL2\n etc (matching cont. time samples)
%
% EmgOut = ['Emg1' 'Emg2' 'Emg3' 'Emg4' ; emg_chan']