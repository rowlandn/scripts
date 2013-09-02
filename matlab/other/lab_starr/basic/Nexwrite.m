% Nexwrite.m
% FIFTH FILE CALLED BY Nexmerge.m
% 
% Completely INDEPENDENT of which of the previous files was toggled on/off
% Saves the current NEX file structures: one "Nexfile." and "Nexvar."
% along with "Nexdata." for each variable/channel.
% NOTE: implement UNIQUE CHANNEL NAME CHECKING !!!

% *******************************************************************
% RE-CALCULATE the OFFSET FIELDS for ALL variables
Nexvar.DataOffset(1,1) = 544 + Nexfile.NumVars*208; % initial offset, just after all Nexfile.NumVars headers
%%% SHOULD "Nexvar.Count(i,:)" be used for offsets/writing ????
for i = 2:Nexfile.NumVars % Cycle through remaining variables
    % WE MUST ONLY WORRY ABOUT 0/1 and 5 for now.. 2/6 soon.. and 3 is EASY.. 4/6 ?????
    switch Nexvar.Type(i-1,1)
        case {0, 1} % set offset if PREVIOUS neuron/event
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i-1,1) + 4*length(Nexdata.TS{i-1,1});
        case 2  % set offset if PREVIOUS interval (ts1, ts2)
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i-1,1) + 8*length(Nexdata.TS{i-1,1});
        case 3  % set offset if PREVIOUS waveform (ts, ad)
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i-1,1) + 4*length(Nexdata.TS{i-1,1});
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i,1) + 2*length(Nexdata.Cont{i-1,1});
        case 4  % set offset if PREVIOUS population vector
            % SKIP THIS FOR NOW
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i-1,1) + 4*length(Nexdata.TS{i-1,1});
            disp(sprintf('NOT SURE HOW TO CALCULATE "population vector" OFFSET !\n'));
        case 5  % set offset if PREVIOUS continuous A/D (ts, ad, index)
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i-1,1) + 4*length(Nexdata.TS{i-1,1});
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i,1) + 4*length(Nexdata.FragIndices{i-1,1});
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i,1) + 2*length(Nexdata.Cont{i-1,1});
        case 6  % set offset if PREVIOUS marker
            Nexvar.DataOffset(i,1) = Nexvar.DataOffset(i-1,1) + 4*length(Nexdata.TS{i-1,1});
            % SKIP THIS FOR NOW
            disp(sprintf('NOT SURE HOW TO CALCULATE "marker" OFFSET !\n'));
    end
end

% ###################################################################
% WRITING the FILE HEADER crap to construct a new NEX file (add "new file")
% ###################################################################
%%%%%%%%%%%%%%%%%%%%%%%%%% Make nice file SAVE dialogue later!
% SHOULD DO:  set the tokens here:  if ~merge_cont, token='', if ~is_event, token='', etc
% could calculate # of events stored by findstr(nexfilename, evt_token)
fnd = findstr(nexfilename, merge_token);
if length(fnd) % move the merge_token to end of the name
    for i=length(fnd):-1:1 % count backwards to avoid screwing up indices
        nexfilename(fnd(i):fnd(i)+length(merge_token)-1) = []; % Delete previous internal merge tokens from filename
    end
end
nexfilename = [nexfilename(1:end-4) merge_token '.nex'];
fid = fopen([nexpath nexfilename], 'w'); % for now, any file with this NAME must be PRE-EXISTING

fwrite(fid, Nexfile.MagicNumber, 'int32');
fwrite(fid, Nexfile.Version, 'int32');
%%% DISPLAY CURRENT CHANNEL NAMES FOR ASSISTANCE EDITING OLD COMMENT
if rename_comment
    disp(sprintf('\n\n\n\nNow preparing to write NEW nex file...\n'));
    disp(sprintf('===================================================================================='));
    disp(sprintf('The current CHANNEL NAMES are:\n'));
    disp(char(Nexvar.Name))
    %%% SPACES SHOULD NOT MATTER, OLD COMMENT IS VERBATIM FROM THE *.nex file
    disp(sprintf('\nThe current NEX FILE COMMENT is:\n\t%s\n', deblank(Nexfile.Comment)))
    choice=input('WOULD YOU LIKE TO CREATE A NEW NEX FILE COMMENT???\n If so, type "y" and press ENTER. If not, just press ENTER\n # > ', 's');
    while ~strcmpi(choice, '')
        if strcmpi(choice, 'y')
            Nexfile.Comment = input('\n\nType NEW nex file comment now.  Feel free to copy and paste from the old one. \n REMEMBER, ONLY 256 CHARACTERS WILL BE WRITTEN!!!\n # > ', 's');
            choice = '';
        else
            choice=input('WOULD YOU LIKE TO CREATE A NEW NEX FILE COMMENT???\n If so, type "y" and press ENTER. If not, just press ENTER\n # > ', 's');
        end
    end
end
commentsize=length(Nexfile.Comment);
if commentsize < 256
    empty = 256 - commentsize;       % How many bytes are empty?
    empty = zeros(1, empty);     % Should be faster than a loop writing '0'
    fwrite(fid, Nexfile.Comment,'uchar');
    fwrite(fid, empty, 'uchar');
else 
    fwrite(fid, Nexfile.Comment(1:256), 'uchar'); % if comment GREATER than 256, truncate it
end
% back to simplicity
fwrite(fid, Nexfile.Frequency, 'float64');
fwrite(fid, Nexfile.Beg, 'int32');
fwrite(fid, Nexfile.End, 'int32');         % FILE HEADER VALUE MAY BE ALTERED !!!
fwrite(fid, Nexfile.NumVars, 'int32');     % FILE HEADER VALUE MAY BE ALTERED !!!
fwrite(fid, Nexfile.NextFileHeader, 'int32');
fwrite(fid, Nexfile.Padding, 'uchar');

% ###################################################################
% WRITING the basic VARIABLE HEADER crap to construct a new NEX file
% ###################################################################

% control this by a loop counter variable (matrix row = channel variable)
disp(sprintf('\n\n\n===================================================================================='));
for i = 1:Nexfile.NumVars
    fwrite(fid, Nexvar.Type(i,:), 'int32');
    fwrite(fid, Nexvar.Version(i,:), 'int32');
    %%% code to WRITE/ALTER variable length channel names: DEBLANK the padded strings!
    varname = Nexvar.Name(i,:);
    if rename_variables
        disp(sprintf('Variable/Channel #%d (type %d) is named:\n\n\t%s\n', i, Nexvar.Type(i,:), deblank(Nexvar.Name(i,:))));
        choice=input('WOULD YOU LIKE TO CHANGE THIS ???\n If so, type "y" and press ENTER. If not, just press ENTER\n # > ', 's');
        while ~strcmpi(choice, '')
            if strcmpi(choice, 'y')
                varname = input('\nType new VARIABLE NAME now.  Feel free to copy and paste. \n REMEMBER, ONLY 64 CHARACTERS WILL BE WRITTEN!!!\n # > ', 's');
                choice = '';
            else
                choice=input('WOULD YOU LIKE TO CHANGE THIS ???\n If so, type "y" and press ENTER. If not, just press ENTER\n # > ', 's');
            end
        end
    end
    varname=deblank(varname);
    varnamesize=length(varname);
    if varnamesize < 64
        empty = 64 - varnamesize;
        empty = zeros(1, empty);
%       Nexvar.Name = char(Nexvar.Name(:,:), varname);
        fwrite(fid, varname, 'uchar');
        fwrite(fid, empty, 'uchar');
    else
        fwrite(fid, varname(1:64), 'uchar');
%         varname = varname(1:64);
%         Nexvar.Name = char(Nexvar.Name(:,:), varname);
    end
    % back to simplicity
    fwrite(fid, Nexvar.DataOffset(i,:), 'int32');
    fwrite(fid, Nexvar.Count(i,:), 'int32');
    fwrite(fid, Nexvar.WireNumber(i,:), 'int32');
    fwrite(fid, Nexvar.UnitNumber(i,:), 'int32');
    fwrite(fid, Nexvar.Gain(i,:), 'int32');
    fwrite(fid, Nexvar.Filter(i,:), 'int32');
    fwrite(fid, Nexvar.XPos(i,:), 'float64');
    fwrite(fid, Nexvar.YPos(i,:), 'float64');
    fwrite(fid, Nexvar.WFrequency(i,:), 'float64');  % ASK ???
    fwrite(fid, Nexvar.ADtoMV(i,:), 'float64');
    fwrite(fid, Nexvar.NPointsWave(i,:), 'int32');   % ASK ???
    fwrite(fid, Nexvar.NMarkers(i,:), 'int32');
    fwrite(fid, Nexvar.MarkerLength(i,:), 'int32');
    fwrite(fid, Nexvar.Padding(i,:), 'uchar');
end

% ###################################################################
% WRITING the ACTUAL CHANNEL DATA to complete the new NEX file
% ###################################################################

% NOW WRITE THE DATA (TS should be mult by Nexfile.Frequency to make it integral time ticks)
% It may be a good idea to just ALTER nex_ts() so it does NOT divide the raw Timestamp Values
for i=1:Nexfile.NumVars
    switch Nexvar.Type(i,:)
        case {0, 1} % write data for neuron/event
            fwrite(fid, Nexfile.Frequency*Nexdata.TS{i,:},'int32');
        case 2  % write data for interval (ts1, ts2)
            fwrite(fid, Nexfile.Frequency*Nexdata.TS{i,:},'int32');
            fwrite(fid, Nexfile.Frequency*Nexdata.TS2{i,:},'int32');
        case 3  % write data for waveform (ts, ad)
            fwrite(fid, Nexfile.Frequency*Nexdata.TS{i,:}, 'int32');
            fwrite(fid, Nexdata.Cont{i,:}, 'int16');  
        case 4  % write data for population vector
            % SKIP THIS FOR NOW
            disp(sprintf('NOT SURE HOW TO WRITE "population vector" !\n'));
        case 5  % write data for continuous A/D (ts, ad, index)
            fwrite(fid, Nexfile.Frequency*Nexdata.TS{i,:}, 'int32');
            fwrite(fid, Nexdata.FragIndices{i,:}, 'int32');
            fwrite(fid, Nexdata.Cont{i,:}/Nexvar.ADtoMV(i,:), 'int16');
        case 6  % write data for marker (ts, names,
            fwrite(fid, Nexfile.Frequency*Nexdata.TS{i,:}, 'int32');
            for j=1:Nexdata.NumFields(i,1)
                fwrite(fid, Nexdata.FieldNames(j,:), 'char'); % WRITE jth FIELD NAME
                for k=1:Nexdata.TS_count(i,1)
                    % WRITE k MARKERS FOR jth FIELD
                    fwrite(fid, Nexdata.Cont{i,1}(k, :, j), 'char'); 
                end
            end
    end
end

fclose(fid);

disp(['successful NEX file re-creation at ' nexpath nexfilename]);