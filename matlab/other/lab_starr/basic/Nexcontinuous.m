% Nexcontinuous.m
% SECOND FILE CALLED BY Nexmerge.m
% 
% USES following parameters, create if necessary!:
if ~exist('automate_fileloads')
    automate_fileloads = 0;
end
if ~exist('automate_variableloads')
    automate_variableloads = 0;
end
%####################################################################

% *******************************************************************
% LOADING *.mat FILES with variables for type 5 CONTINUOUS DATA
% *******************************************************************
loadstring = 'Select MAT file with Continuously Sampled Data, or press CANCEL';
if automate_fileloads
    % AUTOMATED FILE LOADING of matfile CORRESPONDING to nexfile.
    if exist([nexpath basefilename(1:end-4) '.mat']) == 2
        load([nexpath basefilename(1:end-4) '.mat']);
    else
        automate_fileloads = 0;
        loadstring = 'MATFILE NOT FOUND: Select MAT file with Cont Data, or CANCEL';
    end
end
if ~automate_fileloads
    % MANUAL FILE LOADING
    [matfilename, matpathname] = uigetfile('*.mat', loadstring);
    if matfilename == 0 & matpathname == 0   % Allow for exiting out of this script if CANCEL
        n_emg = 0; % a really shitty exit condition!!!
        automate_variableloads = 0; % a really shitty exit condition!!!
    else
        load([matpathname matfilename]);  % GRABS CONTINUOUS DATA: n_emg, emg_chan, time
    end
end

% *******************************************************************
% FETCHING the ACTUAL DATA values for Type 5 Continuous Data
% *******************************************************************
% COMING SOON: PROMPT USER FOR CHANNEL / TIMES(optional) / FREQUENCY
%%% ASSUMES: 1 Nexdata.Cont value / ts since that is the format of our Matlab data!
%%% ASSUMES: "emg_chan" = cont data, "n_emg" = num channels, "time" = timestamps

if automate_variableloads
    if from_abfconv % default:
        % time = time array of 1 X N doubles;
        % n_emg = # of channels;
        % emg_chan = data array of n_emg chans X N doubles;
        % PERHAPS: emg_names = char array of n_emg chans X 64 chars
        if (exist('emg_names') == 1)
            [blah blah2] = size(emg_names);
            if (blah ~= n_emg) | (blah2 ~= 64) % check name array dimensions
                rename_channels=1;
            else
                rename_channels=0;
            end
        else
            rename_channels=1;
        end
        if rename_channels  % Just auto-name variables numerically
            emg_names='';
            for i=1:n_emg
                if i==n_emg % Make "n_emg" ?
                    varname=[num2str(i) '_Accel'];
                else
                    varname=num2str(i);
                end
                varnamesize=length(varname);
                if varnamesize < 64
                    empty = 64 - varnamesize;
                    empty = char(zeros(1, empty));
                    varname = char([varname empty]);
                else
                    varname = varname(1:64);
                end
                emg_names(i,:) = varname
            end
        end
    elseif from_intraop
        % SIMPLY RE-MAP CRAP to make data fit abfconv names:
        [blah n_emg] = size(daq1_deci_data);
        emg_chan = double(daq1_deci_data.'); % transpose to fit our data structure
        time = daq1_deci_times.'; % transpose to fit our data structure
        % RLM: added 2003.10.02, NOT tested.
        emg_names='';
        for i=1:n_emg
            if i==n_emg % Make "n_emg" ?
                varname=[num2str(i) '_Accelerometer'];
            else
                varname=num2str(i);
            end
            varnamesize=length(varname);
            if varnamesize < 64
                empty = 64 - varnamesize;
                empty = char(zeros(1, empty));
                varname = char([varname empty]);
            else
                varname = varname(1:64);
            end
            emg_names(i,:) = varname
        end        
    else
        % ADD FUTURE CASES OR GENERIC FORMAT
    end
else
    % HERE PROMPT FOR VARIABLES IN WORKSPACE
    % GET TIME VAR, GET DATA VAR(s), GET DATA DIMENSIONS, GET CHAN NAMES
    % Map prompted variables to the new names!
    whos % print out workspace variables for selection
end

% *******************************************************************
% LOAD CONTINUOUS DATA into NEX Variable Structure
% *******************************************************************
for i = Nexfile.NumVars+1:Nexfile.NumVars+n_emg  % n_emg (1x1) = number of continuous channels in *.mat file
    Nexdata.TS{i,1} = time;  % time (1xNexvar.Count) = time stamps in DECIMAL SECONDS (x.xxxxx):
    Nexdata.TS_count(i) = length(Nexdata.TS{i,1});
    Nexdata.Freq(i) = round(length(time)/(time(end)-time(1)));
    if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK FILE BOUNDS
        Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
    end
    Nexdata.Cont{i,1} = emg_chan(i-Nexfile.NumVars,:); % LOADS A/D CONTINUOUS VALUES
    % LOAD INDICES INTO TIMESTAMPS (ours are ASSUMED 1 fragment / ts)
    Nexdata.FragIndices{i,1} = 0:length(Nexdata.Cont{i,1})-1
    if Nexdata.FragIndices{i,1}(1,end) > length(Nexdata.Cont{i,1}(1,:))-1
        disp('INCORRECT FRAGMENT INDICES');
    end
    
%  COMPLETE VARIABLE HEADERS FOR NEW TYPE 5 DATA 
    Nexvar.Type(i,1) = 5;        % CONTINUOUS is Type 5
    Nexvar.Version(i,1) = 100;   % seems to be ~constant...
%%% how should cont channels be NAMED?? it's not stored in *.mat file !! % PROMPT FOR A NAME ???
    Nexvar.Name(i,:) = emg_names(i-Nexfile.NumVars,:);
%OLD:    Nexvar.Name = char(Nexvar.Name(:,:), strcat('emg_chan', num2str(i-Nexfile.NumVars))); % use "char()" to add trailing spaces
%%% Nexvar.DataOffset is COMPLICATED (INCREMENTING VALUE IS JUST FOR HEXEDITOR DEBUGGING)
    Nexvar.DataOffset(i,1) = i;  % DUMMY VALUE, we'll calculate far below
%%% Nexvar.Count should be made independent, based upon each separate READ / DATA CHANNEL INPUT
    Nexvar.Count(i,1) = length(Nexdata.TS{i,1});  % typically ALL CONTINUOUS VARIABLES SAME LENGTH
    % next 6 are USELESS (only for neuronal data)
    Nexvar.WireNumber(i,1) = 0;
    Nexvar.UnitNumber(i,1) = 0;
    Nexvar.Gain(i,1) = 0;
    Nexvar.Filter(i,1) = 0;
    Nexvar.XPos(i,1) = 0;
    Nexvar.YPos(i,1) = 0;
    % next 5 are only useful in CONTINUOUS/WAVEFORM/MARKER variables, NOT in neuronal
    if from_abfconv
        Nexvar.WFrequency(i,1) = round(length(time)/(time(end)-time(1)));
    elseif from_intraop
        Nexvar.WFrequency(i,1) = 1000;   % CONTINUOUS DATA SAMPLING RATE.  Get user input, calculate ???
    else
        Nexvar.WFrequency(i,1) = 1000;
    end
% Scale factor to WRITE INTEGRAL VALS !
    maxval = max(abs(Nexdata.Cont{i,1}))
    Nexvar.ADtoMV(i,1) = maxval/32000;  % Autoscale: +/-2^15 = +/-32,768 > 32,000
    Nexvar.NPointsWave(i,1) = length(Nexdata.Cont{i,1}); % CONSTANT IF WAVEFORM (3), total points if CONTINUOUS (5)
    Nexvar.NMarkers(i,1) = 0;        % ask IF MARKER (type 6)
    Nexvar.MarkerLength(i,1) = 0;    % ask IF MARKER (type 6)
    Nexvar.Padding(i,:) = zeros(1, 68)
end

Nexfile.NumVars = Nexfile.NumVars + n_emg;   % set CORRECT number of variables for the file...
if n_emg > 0
    merges = 0;
    fnd = findstr(nexfilename, continuous_token); % is basefile result of previous nex merges?
    if length(fnd)
        for i=length(fnd):-1:1 % count backwards to avoid screwing up indices
            j = str2num(nexfilename(fnd(i)+length(continuous_token)));
            nexfilename(fnd(i):fnd(i)+length(continuous_token)) = []; % Delete previous internal token+# from filename
            merges = max(merges,j); % find maximum value appended to nex_token (# of nex files currently merged)
        end
    end
    merges = num2str(merges+1);
    nexfilename = [nexfilename(1:end-4) continuous_token merges '.nex']; % append token signifying continuous data merged!
    disp('CONTINUOUS DATA SUCCESSFULLY LOADED')
else
    disp('NO CONTINUOUS DATA LOADED')
end