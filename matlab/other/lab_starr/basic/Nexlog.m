% Nexlog.m
% FOURTH FILE CALLED BY Nexmerge.m

if automate_fileloads
    % AUTOMATED FILE LOADING of CORRESPONDING logfile.
    coglog = loadlog([nexpath nexfilename(1:end-4) '.log'])
else
    % MANUAL FILE LOADING
    [cogentfilename, cogentpath] = uigetfile('*.log', 'Select Cogent LOG file with Channel Events');
    coglog = loadlog([cogentpath cogentfilename]);
end

i = Nexfile.NumVars % Cumulative number of variables (including continuous)

% BYPASS THE COGENT HEADER, BYPASS INITIALIZATION CRAP, 
% FIND LOG OFFSET TIME, FIND FIRST 'BEGIN' to set correct log index
for j=1:length(coglog)  % "j" is ALWAYS the index which parses the coglog
    logtime(j,1) = 0;
    logtext{j,1} = '';
    if length(coglog{j}) >= 4
        if strcmp(coglog{j}{4},'COGENT') % FIRST MEANINGFUL LOG ENTRY
            logtime(j,1) = str2num(coglog{j}{1});
            logtext{j,1} = coglog{j}{4};
            log_start_index = j;
        elseif strcmp(coglog{j}{4}, 'DAQ1_TRIGGERED') % FINDS THE LOG OFFSET TIME
            logtime(j,1) = str2num(coglog{j}{1});
            logtext{j,1} = coglog{j}{4};
            logoffset = str2num(coglog{j}{1})/1000 % SUBTRACT this many SECONDS from event/log times
            daq_start_index = j;
        elseif strcmp(coglog{j}{4},'BEGIN') % FIRST TRIAL BEGINS
            logtime(j,1) = str2num(coglog{j}{1});
            logtext{j,1} = coglog{j}{4};
            event_start_index = j;
            break
        end
    end
end

events = {}; % Stores EVENTS map (name to index)
intervals = {}; % Stores INTERVALS map (name to index)
typedstring = []; % Stores NUMERIC value of keypresses during trial
% Extend Variable headers to include intervals/events
Nexvar.Name = char(Nexvar.Name(1:i,:),'','','',''); % append 4 empty names for placeholders?
event_num = size(events,1); % "event_num" SHOULD BE ZERO
interval_num = size(intervals,1); % "interval_num" SHOULD BE ZERO

% CURRENTLY ONLY WORKS WITH ON PRE-DEFINED INTERVALS:
for j=event_start_index : length(coglog) % "j" parses Cogent Log LINE-BY-LINE
    logtime(j,1) = str2num(coglog{j}{1})/1000 - logoffset; % 1st item in each log line is TIME
    logtext{j,1} = coglog{j}{4}; % We "OPERATE" on the 4th item in each log line
    switch logtext{j,1}
        case 'BEGIN' % INTERVALS begin
            typedstring = []; % CLEAR THE KEYBOARD CAPTURE BUFFER
            current_interval_begin = logtime(j); % SAVE TS FOR LATER!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'END'   % INTERVALS end
            % Check if must CREATE BRAND NEW interval name
            if ~strcmp(deblank(Nexvar.Name(i+1,:)),'ALLTRIALS')
                Nexvar.Name(i+1,:) = ['ALLTRIALS' zeros(1,64-length('ALLTRIALS'))]
                Nexvar.Type(i+1,1) = 2;
                Nexdata.TS{i+1,1}(1) = current_interval_begin;
                Nexdata.TS2{i+1,1}(1) = logtime(j);
                interval_num = interval_num + 1;
                intervals{interval_num,1} = i+1;
                intervals{interval_num,2} = 'ALLTRIALS';
            else
                Nexdata.TS{i+1,1}(end+1) = current_interval_begin;
                Nexdata.TS2{i+1,1}(end+1) = logtime(j);
            end
            %SEARCH FOR CREATED INTERVAL NAMES
            found = 0;
            for k=1:interval_num % SEARCH for current event variable
                if strcmp(intervals{k,2}, current_interval)
                    found = found+1;
                    offset_old = intervals{k,1}; %
                    Nexdata.TS{offset_old,1}(end+1) = current_interval_begin; % append start time
                    Nexdata.TS2{offset_old,1}(end+1) = logtime(j); % append end time
                end
            end
            if found ~= 1 % IF SEARCH FAILS, create interval variable
                interval_num = interval_num+1  % create new event
                intervals{interval_num,1} = offset_next;
                intervals{interval_num,2} = current_interval
                Nexvar.Name(offset_next,:) = [current_interval zeros(1,64-length(current_interval))]
                Nexvar.Type(offset_next,1) = 2;
                Nexdata.TS{offset_next,1}(1) = current_interval_begin
                Nexdata.TS2{offset_next,1}(1) = logtime(j) % append end time
                offset_next = offset_next+1
            end
            % RECONSTRUCT TYPED STRINGS, PRINT TO CONSOLE!
            realstring = '';
            for k=1:length(typedstring) % move through numbers
                for L = 1:length(keydata) % search through indices
                    if typedstring(k) == keydata{L,1}
                        realstring = [realstring keydata{L,2}];
                    end
                end
            end
            
            disp(['Keyboard capture: ' realstring]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'go'    % EVENT 1, create variable first time
            if ~strcmp(deblank(Nexvar.Name(i+2,:)),logtext{j,1})
                Nexvar.Name(i+2,:) = [logtext{j,1} zeros(1,64-length(logtext{j,1}))];
                Nexvar.Type(i+2,1) = 1;
                Nexdata.TS{i+2,1}(1) = logtime(j);
                event_num = event_num + 1;
                events{event_num,1} = i+2; % fixed index
                events{event_num,2} = logtext{j,1};
            else
                Nexdata.TS{i+2,1}(end+1) = logtime(j);
            end
        case 'relax' % EVENT 2, create variable first time
            if ~strcmp(deblank(Nexvar.Name(i+3,:)),logtext{j,1})
                Nexvar.Name(i+3,:) = [logtext{j,1} zeros(1,64-length(logtext{j,1}))];
                Nexvar.Type(i+3,1) = 1;
                Nexdata.TS{i+3,1}(1) = logtime(j);
                event_num = event_num + 1;
                events{event_num,1} = i+3; % fixed index
                events{event_num,2} = logtext{j,1};
            else
                Nexdata.TS{i+3,1}(end+1) = logtime(j);
            end
        case 'Key'  % CONVERT KEY# to CHAR
% WORKS!  Now display it nicer
            if strcmp(coglog{j}{6}, 'DOWN')
                typedstring = [typedstring, str2num(coglog{j}{5})]
            end
        case 'FINISHED' % Do nothing now
        case 'DAQ1_STOPPED'
            % Should grab TIME and compare shit for accuracy
            break % CEASE PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        otherwise % NOW create EVENTS "move" "right" "moveright" etc
            if event_num == 0 % FIRST NEW event
                second_name = 0 % Will be auto-switched to "1" below
                event_num = event_num+1 % one new event (THREE PREVIOUS VARIABLES)
                events{1,1} = i+4 % VARIABLE INDEX INTO DATA (Nexdata.TS, Nexvar.Name)
                events{1,2} = logtext{j,1} % Nexvar.Name for event
                offset_next = events{1,1}
                Nexvar.Name(offset_next,:) = [logtext{j,1} zeros(1,64-length(logtext{j,1}))]
                Nexvar.Type(offset_next,1) = 1;
                Nexdata.TS{offset_next,1}(1) = logtime(j)
                offset_next = offset_next + 1
            else % DEFAULT
                found = 0;
                for k=1:event_num % SEARCH for current event variable
                    if strcmp(events{k,2}, logtext{j,1})
                        found = found+1;
                        offset_old = events{k,1}; %
                        Nexdata.TS{offset_old,1}(end+1) = logtime(j); % append time
                    end
                end
                if found ~= 1 % IF NOT FOUND, create it
                    event_num = event_num+1  % create new event
                    events{event_num,1} = offset_next;
                    events{event_num,2} = logtext{j,1}
                    Nexvar.Name(offset_next,:) = [logtext{j,1} zeros(1,64-length(logtext{j,1}))]
                    Nexvar.Type(offset_next,1) = 1;
                    Nexdata.TS{offset_next,1}(1) = logtime(j)
                    offset_next = offset_next+1
                end
            end
            % SET UP INTERVAL NAMES (moveleft, sayright, etc)
            if second_name == 1
                current_interval = [current_interval logtext{j,1}]; % FULL Name of interval
                second_name = 0;
            else
                current_interval = logtext{j,1}; % 1st HALF Name of interval
                second_name = 1;
            end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% DEBLANK TRAILING VARIABLE NAMES (Probably unnecessary in finished code)
lastindex = size(Nexvar.Name,1);
for k=lastindex:-1:1
    while strcmp(deblank(Nexvar.Name(end,:)),'')
        Nexvar.Name(end,:) = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

for i = Nexfile.NumVars+1:length(Nexvar.Name(:,1)) % Previous Num_Vars to current total
    if Nexdata.TS{i,end}(1,end)*Nexfile.Frequency > Nexfile.End % CHECK FILE BOUNDS
        Nexfile.End = Nexdata.TS{i,end}(1,end)*Nexfile.Frequency + 1;
    end
    % LOAD INDICES INTO TIMESTAMPS (ours are ASSUMED 1 fragment / ts)
    
%  COMPLETE VARIABLE HEADERS FOR NEW TYPE 1/2 DATA    
%    Nexvar.Type(i,1) = PRE-ASSIGNED ABOVE
    Nexvar.Version(i,1) = 100;   % seems to be constant...
%    Nexvar.Name = PRE-ASSIGNED ABOVE FROM SOUND FILE
%%% Nexvar.DataOffset is COMPLICATED (INCREMENTING VALUE IS JUST FOR HEXEDITOR DEBUGGING)
    Nexvar.DataOffset(i,1) = i;  % DUMMY VALUE, we'll really calculate it later
%%% Nexvar.Count should be made independent, based upon each separate READ / DATA CHANNEL INPUT
    Nexvar.Count(i,1) = length(Nexdata.TS{i,1});  % typical
    % next 6 USELESS if NOT neuronal data
    Nexvar.WireNumber(i,1) = 0;
    Nexvar.UnitNumber(i,1) = 0;
    Nexvar.Gain(i,1) = 0;
    Nexvar.Filter(i,1) = 0;
    Nexvar.XPos(i,1) = 0;
    Nexvar.YPos(i,1) = 0;
    % next 5 only useful in CONTINUOUS/WAVEFORM/MARKER variables
    Nexvar.WFrequency(i,1) = 0;
    Nexvar.ADtoMV(i,1) = 0;
    Nexvar.NPointsWave(i,1) = 0;
    Nexvar.NMarkers(i,1) = 0;        % ask IF MARKER (type 6)
    Nexvar.MarkerLength(i,1) = 0;    % ask IF MARKER (type 6)
    Nexvar.Padding(i,:) = zeros(1, 68);
end

Nexfile.NumVars = length(Nexvar.Name(:,1));   % set CORRECT number of variables for the file...
disp('COGENT LOGFILE SUCCESSFULLY LOADED')