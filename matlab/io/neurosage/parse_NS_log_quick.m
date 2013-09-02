function parse_NS_log_quick(filename)

%%%%% DEBUG
%filename = 'D:\Lab Sync\Data\Raw\viv06\viv0624a.data';
%%%%%%%%%

%%%%% TO DO
%%% change the way the for loop adds channels to log
%%% write <none> for channels with no spike discrimination parameters
%%%%%%%%%

%% Start reading data file %%%%%%%%%%%%%%%%%%%%%
info_trials = readneurosage(filename,[],[]);

no_trials = info_trials.NumTrials;

info = readneurosage(filename,[1:no_trials],[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create master log cell with all info %%%%%%%%
log{1,1} = 'Seq Name';
log{1,2} = 'Set Name';
log{1,3} = 'Max Rate';
log{1,4} = 'Time Stamp';
log{1,5} = 'Comments';
log{1,6} = 'Acq Channels';
log{1,7} = 'Stim Analog Channels';
log{1,8} = 'Stim Digital Channels';

for i = 1:no_trials
    log{i+1,1} = info.Trials(i).SeqName;
    log{i+1,2} = info.Trials(i).SetName;
    log{i+1,3} = info.Trials(i).MaxRate;
    log{i+1,4} = info.Trials(i).TimeStamp;
    log{i+1,5} = info.Trials(i).Comments;
    no_channels = size(info.Trials(i).AcqChannels,1);
    for j = 1:no_channels
        channels(1,j) = info.Trials(i).AcqChannels(j).Channel;
    end
    channels = sort(channels);
    log{i+1,6} = num2str(channels);
end
assignin('base','log',log)
assignin('base','channels',channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check channels present %%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(log,1)-1
    chan_1(i,1:size(char(log{i+1,6}),2)) = char(log{i+1,6});
end
chan_2 = unique(chan_1,'rows');
assignin('base','chan_2',chan_2)
for i = 1:size(chan_2,1)
    chan_2(i,1:size(chan_2(i,:),2)+2)  = ['  ',chan_2(i,:)];
end
chan_3 = [''];
for i = 1:size(chan_2,1)
    chan_3(1,size(chan_3,2)+1:size(chan_2(i,:),2)+size(chan_3,2)) = chan_2(i,:);
end

chan_mat_1 = str2num(chan_3);
chan_mat_2 = chan_mat_1(:);
chan_mat_2 = sort(chan_mat_2);
chan_mat_3 = unique(chan_mat_2)';
chan_mat_4 = num2str(chan_mat_3);
chan_mat_5 = chan_mat_3';
assignin('base','chan_mat_1',chan_mat_1)
assignin('base','chan_mat_2',chan_mat_2)
assignin('base','chan_mat_3',chan_mat_3)
assignin('base','chan_mat_4',chan_mat_4)
assignin('base','chan_mat_5',chan_mat_5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Write string variable %%%%%%%%%%%%%%%%%
cell_size = cellfun('size', log, 2);
for i = 1:size(log,2)-1
    cell_sum(i,1) = cell_size(i+1,1) + cell_size(i+1,2) + cell_size(i+1,6);
end
cell_max = max(cell_sum) + 5;

find_filename_01 = strfind(filename,'\');
find_filename_02 = strfind(filename,'.data');
find_filename_03 = filename(find_filename_01(end)+1:find_filename_02-1);
size_filename_03 = size(find_filename_03,2);
asterisks = '******************************';
size_asterisks = size(asterisks,2);
top_line_num_spaces_needed = cell_max - [size_filename_03+size_asterisks];
top_line = [asterisks,blanks(1),find_filename_03,blanks(top_line_num_spaces_needed)];
parse_string(1,:) = top_line;
size_parse_string = size(parse_string,2);

for i = 1:size(log,1)
    if i == size(log,1)
        str_length = size([num2str(i),blanks(2),char(log{i,1}),blanks(2),char(log{i,6})],2);
        num_spaces_needed = size_parse_string - str_length;
        parse_string(i+1,:) = [num2str(i),blanks(2),char(log{i,1}),blanks(2),char(log{i,6}),blanks(num_spaces_needed)];
    else
        str_length = size([num2str(i),blanks(2),char(log{i+1,1}),blanks(2),char(log{i+1,6})],2);
        num_spaces_needed = size_parse_string - str_length;
        parse_string(i+1,:) = [num2str(i),blanks(2),char(log{i+1,1}),blanks(2),char(log{i+1,6}),blanks(num_spaces_needed)];
    end
end
size_parse_string_rows = size(parse_string,1);
end_line_num_spaces_needed = size_parse_string - size_asterisks;
parse_string(size_parse_string_rows,:) = [asterisks,blanks(end_line_num_spaces_needed)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Display in editable text box %%%%%%%%%%% 
figure
set(gcf,'Position',[150 40 1000 700],'MenuBar','none','Name','NeuroSage ''.data'' file parser');
Parse_NS_box = uicontrol(gcf,'Style','edit','Position',[100 100 800 500],...
    'String',parse_string,'Max',10,'Min',1,'HorizontalAlignment','left');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%