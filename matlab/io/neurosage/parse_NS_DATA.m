function parse_NS_log(filename)

%filename = '/Raw/viv06/viv0603e.data'

%%%%% TO DO
%%% change the way the for loop adds channels to log
%%% write <none> for channels with no spike discrimination parameters
%%% sort comment trial #'s
%%%%%%%%%

%% Check whether file already exists %%%%%%%%%%%
find_filename_01 = strfind(filename,'\');
find_filename_02 = strfind(filename,'.data');
find_filename_03 = filename(find_filename_01(end)+1:find_filename_02-1);
txt_file = [find_filename_03,'.txt'];
cd('D:\Lab Sync\Data\Analyzed\Level_I\viv06');
LS = ls;

for i = 1:size(LS,1)
    if isempty(findstr(LS(i,:),txt_file)) == 0 
        warndlg(['The file ',txt_file,...
            ' already exists.  Please manually remove the file from the directory and re-execute the function.'],...
            'File already exists!')
        return
    else
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

channel_quest_box = 'Channel(s) Detected';
channel_answer_1 = questdlg(['The following channels were detected: ',chan_mat_4,'. Were the sources of any of the channels changed during acquisition?'],...
    channel_quest_box,'No','Yes','No');
if channel_answer_1 == 'No'
    for i = 1:size(chan_mat_3,2)
        Label_channels_parameters_str = ['Label_channels_parameters{1,',num2str(i),'} = [''Channel ',num2str(chan_mat_5(i)),'''];'];
        eval(Label_channels_parameters_str)
    end
    assignin('base','Label_channels_parameters',Label_channels_parameters)
    Label_channels_dlgbox_title = 'Label Channels';
    Label_channels_userinpt = inputdlg(Label_channels_parameters,Label_channels_dlgbox_title);
    %assignin('base','Label_channels_userinpt',Label_channels_userinpt)
end

for i = 1:size(Label_channels_userinpt,1)
    if isempty(Label_channels_userinpt{i})
        Label_channels_userinpt{i} = '<none>';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Collect online comments %%%
[com_com,com_last,com_ind] = unique(log(:,5));
com_com(1,:) = [];
com_match_01 = strmatch('Comments',com_com,'exact');
com_com(com_match_01,:) = [];

diff_com_ind = diff(com_ind);
find_diff_com_ind = find(diff_com_ind);
find_diff_com_ind = find_diff_com_ind+1;
char_log_com = char(log{:,5});

if isempty(com_com) == 1
else
    for i = 1:size(com_com,1)
        com_match = strmatch(com_com{i},char_log_com,'exact');
        min_com_match = min(com_match);
        comments{i,1} = ['Trial ',num2str(min_com_match),': "',com_com{i},'"'];
    end
end

% assignin('base','com_com',com_com)
% assignin('base','com_last',com_last)
% assignin('base','com_ind',com_ind)
% assignin('base','diff_com_ind',diff_com_ind)
% assignin('base','find_diff_com_ind',find_diff_com_ind)
% assignin('base','comments',comments)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Begin series of queries %%%%%%%%%%%%%%%%%%%%

%% Ask if continuation of another recording %%%
same_neuron_quest_box = 'Continuation';
same_neuron_answer_1 = questdlg(['Is this file a continuation of a log file from a previously recorded neuron?'],...
    same_neuron_quest_box,'No','Yes','No');
if strcmp(same_neuron_answer_1,'No') == 1
else
    same_neuron_answer_2_title = 'Continuation';
    same_neuron_answer_2_parameters = {'Which neuron is this a continuation of?'};
    same_neuron_answer_2 = inputdlg(same_neuron_answer_2_parameters,same_neuron_answer_2_title);
    %assignin('base','same_neuron_answer_2',same_neuron_answer_2)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ask if user wants to add offline comments %%%
    add_offline_quest_box = 'Offline Comments';
    add_offline_com_answer = questdlg(['Would you like to add offline comments?'],add_offline_quest_box,'Yes','No','Yes');
    if strcmp(add_offline_com_answer,'Yes') == 1
        offline_com_parameters = {'Add offline comments below:'};
        offline_com_title = 'Offline comments';
        offline_com_userinpt = inputdlg(offline_com_parameters,offline_com_title,6);
        assignin('base','offline_com_userinpt',offline_com_userinpt)
    else
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start writing text file %%%%%%%%%%%%%%%%%%%%%%
if isunix == 1
    % % Write .com file
    % Parse filename from user input
    find_filename_01 = strfind(filename,'/');
    find_filename_02 = strfind(filename,'.data');
    find_filename_03 = filename(find_filename_01(end)+1:find_filename_02-1)
    com_file = [find_filename_03,'_com'];
    ps_file = [find_filename_03,'.ps'];
    pdf_file = [find_filename_03,'.pdf'];

    % cd to directory where file will be saved
    cd_str = ['cd D:\Data\Raw\viv06\',...
        filename(find_filename_01(end-1)+1:find_filename_01(end-1)+6)];
    assignin('base','cd_str',cd_str)
    eval(cd_str)


    %% Begin writing file %%%%%%
    fid = fopen(com_file,'w');
    text_asterisks = '******************************';
    text_filename = find_filename_03;
    fprintf(fid,'%s \s\s',text_asterisks);
    fprintf(fid,'%s \n',text_filename);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:no_trials
        fprintf(fid,'%s \t',num2str(i));
        fprintf(fid,'%s \t',char(log{i+1,1}));
        fprintf(fid,'%s \n',char(log{i+1,6}));
    end


    fclose(fid)
    ps_str = ['!a2ps -B -R --columns=1 ', com_file,' -o ',ps_file];
    eval(ps_str)
    ps_str = ['!ps2pdfwr ',ps_file,' ',pdf_file];
    eval(ps_str)
    acroread_str = ['!acroread ',pdf_file,'&'];
    eval(acroread_str)
    ps_str_2 = ['!rm ',ps_file];
    eval(ps_str_2)
    com_str = ['!rm ',com_file];
    eval(com_str)
    % com_str_2 = ['!rm com_file2'];
    % eval(com_str_2)
else
    find_filename_01 = strfind(filename,'\');
    find_filename_02 = strfind(filename,'.data');
    find_filename_03 = filename(find_filename_01(end)+1:find_filename_02-1);
    txt_file = [find_filename_03,'.txt'];
    
    cd(['D:\Lab Sync\Data\Analyzed\Level_I\',filename(find_filename_01(end-1)+1:find_filename_01(end-1)+6)])
    % Begin writing file
    fid = fopen(txt_file,'w');
    text_asterisks = '******************************';
    text_filename = find_filename_03;
    fprintf(fid,'%s \s\s',text_asterisks);
    fprintf(fid,'%s \n\n',text_filename);
    
    %% Add online comments %%%%
    if isempty(com_com)
        fprintf(fid,'%s \n\n','**** ONLINE USER COMMENTS ****');
        fprintf(fid,'%s \n\n',text_asterisks);
    else    
        fprintf(fid,'%s \n\n','**** ONLINE USER COMMENTS ***');
        for i = 1:size(comments,1)
            fprintf(fid,'%s \n\n',char(comments(i,:)));
        end
        fprintf(fid,'%s \n\n',text_asterisks);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Add offline comments %%%%
    if strcmp(add_offline_com_answer,'Yes') == 1
        fprintf(fid,'%s \n\n','**** OFFLINE USER COMMENTS ***');
        for i = 1:size(char(offline_com_userinpt),1)
            fprintf(fid,'%s \n\n',char(offline_com_userinpt{1}(i,:)));
        end
        fprintf(fid,'%s \n\n',text_asterisks);
    else 
        fprintf(fid,'%s \n\n','**** OFFLINE USER COMMENTS ***');
        fprintf(fid,'%s \n\n',text_asterisks);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Add channels %%%%%%%%%%%%
    fprintf(fid,'%s \n\n','**** CHANNEL SOURCES *********'); 
    for i = 1:size(Label_channels_userinpt,1)
        if i == 1
            fprintf(fid,'%s \n',[num2str(chan_mat_5(i)),' - ',char(Label_channels_userinpt(i))]);
        elseif i == size(Label_channels_userinpt,1)
            fprintf(fid,'%s \n\n',[num2str(chan_mat_5(i)),' - ',char(Label_channels_userinpt(i))]);
        else
            fprintf(fid,'%s \n',[num2str(chan_mat_5(i)),' - ',char(Label_channels_userinpt(i))]);
        end
    end
    fprintf(fid,'%s \n\n',text_asterisks);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Add spike discrimination parameters %%
    spike_discrim_quest_box = 'Spike Discrimination Parameters';
    spike_discrim_answer_1 = questdlg(['Would you like to save spike discrimination parameters?'],...
    spike_discrim_quest_box,'Yes','No','Yes');
    if strcmp(spike_discrim_answer_1,'Yes') == 1
        Label_channels_spike_discrim_params_dlgbox_title = 'Label Spike Discrimination Parameters';
        Label_channels_spike_discrim_params_userinpt = inputdlg(Label_channels_userinpt,Label_channels_spike_discrim_params_dlgbox_title);
        assignin('base','Label_channels_spike_discrim_params_userinpt',Label_channels_spike_discrim_params_userinpt)
        for i = 1:size(Label_channels_userinpt,1)
            if isempty(Label_channels_spike_discrim_params_userinpt{i,1})
                Label_channels_and_spike_discrim_params_str = ['Label_channels_and_spike_discrim_params{',num2str(i),',1} = [num2str(chan_mat_5(',num2str(i),'))  '' - '' char(Label_channels_userinpt(',num2str(i),'))];'];
                eval(Label_channels_and_spike_discrim_params_str)
            else
                Label_channels_and_spike_discrim_params_str = ['Label_channels_and_spike_discrim_params{',num2str(i),',1} = [num2str(chan_mat_5(',num2str(i),')) '' - '' char(Label_channels_userinpt(',num2str(i),')),'':'',char(Label_channels_spike_discrim_params_userinpt(',num2str(i),'))];'];
                eval(Label_channels_and_spike_discrim_params_str)
            end
        end
        
        fprintf(fid,'%s \n\n','**** SPIKE PARAMETERS ********'); 
        for i = 1:size(Label_channels_and_spike_discrim_params,1)
            if i == 1
                fprintf(fid,'%s \n',char(Label_channels_and_spike_discrim_params(i,:)));
            elseif i == size(Label_channels_and_spike_discrim_params,1)
                fprintf(fid,'%s \n\n',char(Label_channels_and_spike_discrim_params(i,:)));
            else
                fprintf(fid,'%s \n',char(Label_channels_and_spike_discrim_params(i,:)));
            end
        end
        fprintf(fid,'%s \n\n',text_asterisks);    
        %assignin('base','Label_channels_and_spike_discrim_params',Label_channels_and_spike_discrim_params)
    else
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Print trials %%%%%%%%%%%
fprintf(fid,'%s \n\n','**** TRIALS ******************'); 
for i = 1:no_trials
    fprintf(fid,'%s \t',num2str(i));
    fprintf(fid,'%s \t',char(log{i+1,1}));
    fprintf(fid,'%s \n',char(log{i+1,6}));
end
fprintf(fid,'%s \n\n',text_asterisks);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignin('base','fid',fid)
    

%% Close and save file %%%%
%fclose(fid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     ps_str = ['!a2ps -B -R --columns=1 ', com_file,' -o ',ps_file];
%     eval(ps_str)
%     ps_str = ['!ps2pdfwr ',ps_file,' ',pdf_file];
%     eval(ps_str)
%     acroread_str = ['!acroread ',pdf_file,'&'];
%     eval(acroread_str)
%     ps_str_2 = ['!rm ',ps_file];
%     eval(ps_str_2)
%     com_str = ['!rm ',com_file];
%     eval(com_str)
    % com_str_2 = ['!rm com_file2'];
    % eval(com_str_2)
end
    