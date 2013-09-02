
% ################################################################################
% ################################################################################
% RLM: UGLY "merge" function for TWO nex files
% RLM: Still needs to check for UNIQUE NAMES!
% RLM: LATER add "n-files" so we can merge a zillion files if entered in the gui.
if merge_nexfiles
  while 1 % start a loop block, just so we can "break out" by choosing "cancel"
    [nexfilename2, nexpath2] = uigetfile('*.nex', 'Select SECOND NEX FILE to merge INTO the first, OR CANCEL');
    if nexfilename2 == 0 & nexpath2 == 0
        merge_nexfiles = 0;
        is_event = 0; % just to make sure not set already
        break;
    end
    nid2=fopen([nexpath2 nexfilename2],'r');     % Open selected file READ ONLY
    if (nid2 == -1)
        error('Cannot open SECOND NEX file!');
        return
    end
    
    temp_nexname2 = strrep(nexfilename2,event_is_neuron_token,event_was_neuron_token); % convert filename
    fnd = findstr(temp_nexname2,event_was_neuron_token); % look for token
    if length(fnd)
    % IF MERGING "event stored as neuron" NEXFILE (usually accelerometer deflection from Offline Sorter)
        is_event = 1;   % NEW merged Nexfile will have this file's data as "event" NOT "neuron"
            isi_calculations = 0; % Avoid doing any other Spontmerge tasks if importing EVTS
            write_nex = 1; % Duh, of course we will be writing the file.
        % modify the above nexfilename, appending it's "event_was_neuron_token" with # (ex: _b2)
        incr = temp_nexname2(fnd(1)+length(event_was_neuron_token)) % fetch this for naming event variables below
        nexfilename = [nexfilename(1:end-4) temp_nexname2(fnd(1):fnd(1)+length(event_was_neuron_token)) '.nex'];
    else
    % These files are simple merging of 2 regular NEX files!
        is_event = 0;
        merges = 0;  % for FIRST TIME merging of a nex file into the basefile
        fnd = findstr(nexfilename, nex_token); % is basefile result of previous nex merges?
        if length(fnd)
            for i=length(fnd):-1:1 % count backwards to avoid screwing up indices
                j = str2num(nexfilename(fnd(i)+length(nex_token)));
                nexfilename(fnd(i):fnd(i)+length(nex_token)) = []; % Delete previous internal token+# from filename
                merges = max(merges,j); % find maximum value appended to nex_token (# of nex files currently merged)
            end
        end
        merges = num2str(merges+1); % increment to account for the latest merging!
        nexfilename = [nexfilename(1:end-4) nex_token merges '.nex'];
    end

    % *******************************************************************
    % READING/LOADING *.nex file FILE HEADER fields (add more GUI / Error Checking)
    % *******************************************************************
    % NOTE: after READING the file at X precision -> ML converts to DOUBLE or CHAR matrix  
    Nexfile.MagicNumber2 = fread(nid2, 1, 'int32');
    if Nexfile.MagicNumber2 ~= 827868494
        error('This file is NOT a NEX file!');
        return
    end
    Nexfile.Version2 = fread(nid2, 1, 'int32');
    if Nexfile.Version2 ~= Nexfile.Version
        warning('The 2 NEX files have different NEX file versions');
    end
    Nexfile.Comment2(1,:) = fread(nid2, 256, 'uchar=>char'); % READ as 256 unsigned chars -> char
    Nexfile.Frequency2 = fread(nid2, 1, 'float64');
    if ~is_event & (Nexfile.Frequency2 ~= Nexfile.Frequency) % NOTE: ought to have the same timebase for somethings!
        warning('The 2 NEX files have different samplerates!');
%         error('The 2 NEX files have different samplerates!');
%         return
    end
    Nexfile.Beg2 = fread(nid2, 1, 'int32');
    if Nexfile.Beg2 < Nexfile.Beg
        Nexfile.Beg = Nexfile.Beg2;
    end
    Nexfile.End2 = fread(nid2, 1, 'int32');
    if Nexfile.End2 ~= Nexfile.End
        warning('The 2 NEX files have different end times');
        if Nexfile.End2 > Nexfile.End
            Nexfile.End = Nexfile.End2;
        end
    end
    Nexfile.NumVars2 = fread(nid2, 1, 'int32');
    Nexfile.NextFileHeader2 = fread(nid2, 1, 'int32');
    Nexfile.Padding2 = fread(nid2, [1, 256], 'uchar=>char');
    fclose(nid2); % EXIT FILE READING (after reading file header)
    
    % USE NEX AUTHOR'S MATLAB FUNCTION TO GRAB: #vars, varnames, vartypes:
    % THIS ACTUALLY GRABS ALL VARIABLES!! (not just type 0/1)
    [nex_nvar2 nex_names2 nex_types2] = nex_info([nexpath2 nexfilename2]); % grab NUM vars, NAMES vars, TYPES vars
    %%%% IF THE NEX FILE ALREADY HAS NON-NEURONAL VARIABLES, ALERT THE USER !!
    non_ts_indices2 = find(nex_types2 > 1);
    if length(non_ts_indices2 > 0)
        disp(sprintf('Includes non-neuron/non-event variables of type: %d !\n', nex_types2(non_ts_indices2)));
    end
    %clear Nexfile.*2; % Clears the 2nd file header vars

    % *******************************************************************
    % LOADING *.nex file VARIABLE HEADER field variables
    % *******************************************************************
    
    % NOW BEGIN VECTORIZING (or CELLS) NEX VARIABLE HEADERS:
    % Keeps "reopening" the file because of some strange loopspace error...
    for i = Nexfile.NumVars+1:Nexfile.NumVars+nex_nvar2
        readoffset = 544 + (i-(Nexfile.NumVars+1))*208;
        nid=fopen([nexpath2 nexfilename2],'r');
        fseek(nid2, readoffset, 'bof');
        % READ EXACT COPY OF EACH VARIABLE HEADER
        Nexvar.Type(i,:) = fread(nid2, 1, 'int32'); % move filepointer forwards
        if is_event & Nexvar.Type(i,:) == 0 % IF ACCELEROMETER DEFLECTION or other EVENT wrongly stored as neuron
            Nexvar.Type(i,:) = 1; % alter type to be EVENT
        end
        Nexvar.Version(i,:) = fread(nid2, 1, 'int32'); % currently seems to be a constant "100"
        Nexvar.Name(i,:) = fread(nid2, 64, 'uchar=>char'); % IF EVENT, rename below
        %%%% Nexvar.DataOffset CALCULATED LATER (INCREMENTED VALUE JUST FOR HEXEDITOR DEBUGGING)
        Nexvar.DataOffset(i,:) = fread(nid2, 1, 'int32'); %  [REALLY: index of 1st associated data point in file]
        Nexvar.Count(i,:) = fread(nid2, 1, 'int32'); % CHECK THIS FOR ACCURACY LATER
        % BELOW ALL IS FINE FOR NOW:
        % next 6 are usually useless, unless reading in *.nex file with 3-D neuron data
        Nexvar.WireNumber(i,:) = fread(nid2, 1, 'int32');
        Nexvar.UnitNumber(i,:) = fread(nid2, 1, 'int32');
        Nexvar.Gain(i,:) = fread(nid2, 1, 'int32');
        Nexvar.Filter(i,:) = fread(nid2, 1, 'int32');
        Nexvar.XPos(i,:) = fread(nid2, 1, 'float64');
        Nexvar.YPos(i,:) = fread(nid2, 1, 'float64');
        % next 5 are only useful in CONTINUOUS/WAVEFORM/MARKER variables, NOT in neuronal/event
        Nexvar.WFrequency(i,:) = fread(nid2, 1, 'float64');
        Nexvar.ADtoMV(i,:) = fread(nid2, 1, 'float64');
        Nexvar.NPointsWave(i,:) = fread(nid2, 1, 'int32');
        % the next 76 bytes are ALWAYS BLANK (make '0') UNLESS using MARKERS
        Nexvar.NMarkers(i,:) = fread(nid2, 1, 'int32');
        Nexvar.MarkerLength(i,:) = fread(nid2, 1, 'int32');
        Nexvar.Padding(i,:) = fread(nid2, [1, 68], 'uchar=>char');
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
                [Nexdata.TS_count(i,1), Nexdata.TS{i,1}] = nex_ts([nexpath2 nexfilename2], Nexvar.Name(i,:));   % should create Nx64 char arrays
                if is_event & Nexvar.Type(i,:) == 1 % if event in special file, change the name...
                    temp = [def_event_name incr num2str(i-Nexfile.NumVars)]; % name + number
                        % As set in base script (ex: spontmerge.m)
                    Nexvar.Name(i,1:length(temp)) = temp;
                end
                max(Nexdata.TS{i,1})
                Nexdata.Freq(i,1) = Nexfile.Frequency2;
                if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                    Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
                end
                if Nexdata.TS_count(i,1) ~= Nexvar.Count(i,1)
                    disp(sprintf('HEADER %d timestamps not equal to ACQUIRED %d timestamps!\n', Nexvar.Count(i,1), Nexdata.TS_count(i,1)));
                end
            case 2  % READ INTERVAL: beginning & end TIMESTAMP ARRAYS
                [Nexdata.TS_count(i,1), Nexdata.TS{i,1}, Nexdata.TS2{i,1}] = nex_int([nexpath2 nexfilename2], Nexvar.Name(i,:));
                Nexdata.Freq(i,1) = Nexfile.Frequency;
                if Nexdata.TS2{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                    Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
                end
            case 3  % READ WAVEFORM: (ts, ad)
                [Nexdata.Freq(i,1), Nexdata.TS_count(i,1), Nexdata.TS{i,1}, Nexdata.FragSizes{i,1}, Nexdata.Cont{i,1}] ...
                    = nex_wf([nexpath2 nexfilename2], Nexvar.Name(i,:));
                % Nexdata.Cont{i,1} IS NF AD/timestamps x N timestamps/waveforms
                if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                    Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
                end
            case 4  % READ POPULATION VECTOR
                % SKIP THIS FOR NOW
            case 5  % READ CONTINUOUS A/D: ts, 
                [Nexdata.Freq(i,1), Nexdata.TS_count(i,1), Nexdata.TS{i,1}, Nexdata.FragSizes{i,1}, Nexdata.Cont{i,1}] = nex_cont([nexpath2 nexfilename2], Nexvar.Name(i,:));
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
                    = nex_marker([nexpath2 nexfilename2], Nexvar.Name(i,:));
                % FOR MARKERS: Nexdata.Cont{i,1} contains N x NL x NM array
                % n = num markers x chars/field x fields/marker
                % nm = fields/marker, nl = chars/field
                if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK TIME BOUNDS
                    Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
                end
        end
    end
    Nexfile.NumVars = Nexfile.NumVars + nex_nvar2
    disp('SECOND NEX FILE SUCCESSFULLY INTEGRATED');
    break; % get out of block
  end
end