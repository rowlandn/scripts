% Rory L. Michaelis
% 2003.10.09 2003.07.02 2002.12.08 (2002.08.26) (2002.08.19 began)
% 
% SCRIPT for joining pre-existing NEX Data File (*.nex) ...
% (COMING SOON, optionally do this in BATCH MODE)
% ...with optionally another *.nex file
% ...with optionally CONTINUOUS Sampled Data loaded into the Matlab workspace
%   (Cont Sampled Data typically comes from the "abfconv.m" script by Rob Turner)
% ...with optionally CogentLog EVENTS & INTERVALS (hard-coded structure)
% ...with optionally a CogentLog keystroke MARKERS
%
% NOTE: if using HEX EDITOR: x86 architecture uses LITTLE ENDIAN (01 00 = 1, 00 01 = 256)
% -------------------------------------------------------------------
clear all

global current_path; % Path to THISFILE !
global nexpath; % Path to NEXFILE
global basefilename; % store the original nexfile name, used for loading *.mat, *.log file, WILL NOT CHANGE
global nexfilename; % CURRENT NAME of CONSTRUCTED NEXFILE (is modified during the merging process)
global NUMFILES; % total number of NEX files involved in script
% Next one unnecessary ???
global nid; % NEX FILE, fid

% FILE NAMING PARAMETERS
merge_token = '_m'; % RLM: final filename appendage, denotes ANY merged NEX file (one NOT created by a Plexon Inc app)
nex_token = '_n'; % RLM: used for denoting a file resulting from merging 2 NORMAL NEX files, NOT "event_is_neuron" nex files
continuous_token = '_c'; % RLM: appended to filename after merging of any continuous data
event_is_neuron_token = '_a'; % RLM: denotes NEX files incorporating events wrongly stored as neurons!
    % typically accelerometer deflection events that were sorted by OFFLINESORTER
event_was_neuron_token = '_b'; % RLM: fudge added to differentiate file AFTER the neuron-to-event conversion!
def_event_name = 'AccelEvt'; % if merging above type file, uses this default NEX event variable name

% PARAMETER TOGGLES:
if isempty(findobj('Tag','NexGUI')) % If NOT USING GUI
    merge_nexfiles = 1; % Optionally integrate a SECOND Nexfile's data INTO the first Nexfile
        merge_number = 2; % CURRENTLY UNUSED (maybe support crazy mega-merge function later)
    automate_fileloads = 1;  % Assumes *.MAT and *.LOG files use same name as the *.NEX file
    isi_calculations = 0; % Performs ISI calculations on all spike channels, then writes textfile
        isi_batching = 0; % Performs ISI calculations on ALL NEX FILES in directory, writes !!!
        isi_use_threshold = 1; % drop ISI's below a certain threshold
            isithreshold = 1; % IN MILLISECONDS: ISI threshold
            binwidth = 5; % IN MILLISECONDS: width of ISI bins for histogram etc
            normfactor = 10; % IN MILLISECONDS: desired mean ISI when normalized spike train!
    merge_continuous = 0; % Load the *.MAT file containing continuous channels ???
        automate_variableloads = 1; % Assumes CONTINUOUS variables in workspace use PREDEFINED NAMES:
            % from abfconv OR from intraop:  'n_emg', 'emg_chan', 'time', NEW: 'emg_names'
            from_abfconv = 1; % Loading continuous ABF data from RST's abfconv routine
            from_intraop = 0; % Loading continuous NI-DAQ data from IntraOperativeRecording routine
    merge_cogentlog = 0; % Load the Cogent2000 logfile containing events & intervals
        merge_cogentlog_keystrokes = 0; % Store the actual KEYSTROKES during a trial as markers
    write_nex = 1;
        rename_comment = 1; % Prompt user to rename NEX file comment
        rename_variables = 1; % Prompt user to rename all NEX variables before writing new NEX file.
% ############################ BELOW SHOULD NOT BE MODIFIED !!! ################################
else
    % Fetch that shit from the GUI
    disp('Grabbing parameters from NEX gui... ');
    merge_nexfiles = get(findobj('Tag','Chk_MergeNexfiles'),'Value')
        merge_number = str2num(get(findobj('Tag','Edt_NumNexMerged'),'String'))
    automate_fileloads = get(findobj('Tag','Chk_AutomateFileloads'),'Value')
    isi_calculations = get(findobj('Tag','Chk_ISICalculations'),'Value')
        isi_batching = get(findobj('Tag','Chk_ISIBatching'),'Value')
        isi_use_threshold = get(findobj('Tag','Chk_ISIUseThreshold'),'Value')
            isithreshold = str2num(get(findobj('Tag','Edt_ISIMinimumThreshold'),'String'))
    merge_continuous = get(findobj('Tag','Chk_MergeContinuous'),'Value')
        automate_variableloads = get(findobj('Tag','Chk_AutomateVariableLoads'),'Value')
            % from abfconv OR from intraop:  'n_emg', 'emg_chan', 'time', NEW: 'emg_names'
            from_abfconv = 1
            from_intraop = 0
    merge_cogentlog = get(findobj('Tag','Chk_MergeCogentLog'),'Value')
        merge_cogentlog_keystrokes = get(findobj('Tag','Chk_MergeCogentLogKeystrokes'),'Value')
    write_nex = get(findobj('Tag','Chk_WriteNex'),'Value')
        rename_comment = get(findobj('Tag','Chk_RenameComment'),'Value')
        rename_variables = get(findobj('Tag','Chk_RenameVariables'),'Value')
    % variable = get(findobj('Tag','   '),'Value')
    set(findobj('Tag','NexGUI'),'HandleVisibility','Callback');
    disp(' ...Done!');
end

%################################################################################
% -------------------------------------------------------------------
% DOING ACTION TIME:
% Call different scripts based upon the toggles set above
% -------------------------------------------------------------------
%################################################################################

nexinit;        % CALLS M FIILE #0, REQUIRED, sets up things, finds NEX file(s).

if isi_batching & isi_calculations
    for n=1:NUMFILES    % Batch process Nexfile ISI data without merge prompts
        nexload;        % CALLS M FILE #1, REQUIRED, loads pre-existing *.NEX FILE (1)
        nexisi;         % CALLS new m file, OPTIONAL, calculates BATCH ISI stats and writes textfile
    end
    return  % EARLY EXIT, WE JUST WANTED TO SUMMARIZE ALL ISI DATA TO TEXT FILE!
% elseif nex_batching
elseif merge_nexfiles
    nexload;         % CALLS M FILE #1, REQUIRED, loads pre-existing *.NEX FILE (1)
    nexload2;        % CALLS M FILE #1.5, OPTIONAL, attempt to merge 2nd *.NEX FILE with first.
else
    nexload;         % CALLS M FILE #1, REQUIRED, loads pre-existing *.NEX FILE (1)
end
% ######## BELOW THIS LINE SCRIPTS ALWAYS RUN IN NON-BATCH MODE: ##############
% HERE IS WHERE TO PLACE ANALYSIS OF NEX DATA SCRIPTS (single file, non-batched)
if isi_calculations
    nexisi        % CALLS new m file, OPTIONAL, calculates file ISI stats and writes textfile
end
% END ANALYSIS OF NEX DATA SCRIPTS

if merge_continuous
    nexcontinuous;   % CALLS M FILE #2, OPTIONAL, loads continuous channels from *.MAT file or workspace
end

if merge_cogentlog
    if merge_cogentlog_keystrokes
        cogentkeymap    % CALLS M FILE #3, OPTIONAL, loads table of corresponding key numbers & values
    end
    nexlog          % CALLS M FILE #4, OPTIONAL, loads INTERVALS & EVENTS from Cogent *.log file
end

if write_nex
    nexwrite;        % CALLS M FILE #5, OPTIONAL, WRITES new *.NEX FILE with integrated data
end

% ########################################################################################
%           JUST SOME NOTES TO MYSELF:
% *******************************************************************
% NOTE:  when WRITING variables to NEX file, use these precisions:
% 1byte -> 'uchar' 4bytes -> 'int32' 8bytes -> 'float64'
% save (simple), fprintf (versatile), dlmwrite, diary (script copies screen output)

% NeuroExplorer App: FORMAT FOR IMPORTING DATA:  
%   NEURON MULTI-COLUMN-TS Text :  NAMEa\t NAMEb\n TS1a\t TS1b\n TS2a\tTS2b\n
%   OR  LIST Text: NAMEa\t TS1a \n NAMEb\t TS1b \n NAMEb\t  TS2b \n NAMEb\t TS3b \n
%   CONTINUOUS Text: NAME\n VAL1\n VAL2\n etc (matching cont. time samples)
