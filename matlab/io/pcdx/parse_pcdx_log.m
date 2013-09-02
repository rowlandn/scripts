% %%%% TO DO
% 1) Replace channel and channel labels after entering label comments (like spike discrimination parameters)
% with for loop for sed commands
% 
% 2) Make script recognize transitions.  Look at
% /F/Matlab/scripts/tom_log_files/dcn0544b.log for an example.
%     you will just have to go by the pulse delay, width and voltage to know what's going on.  Can't
%     rely on the string describing the transition.  (cip means current injection pulse).


function parse_pcdx_log(filename)
% clear
% !rm /F/Data/Analyzed/Level_I/viv05/viv0518c.pdf
% filename = '/F/Data/Raw/viv05/viv0518c.log';
% This function parses pcdx log files.  Just give the
% '.log' path and filename
% 
% parse_pcdx_log('filename')       <nrowlan@emory.edu>



%clear all
%filename = '/F/Data/Raw/viv05/viv0513h.log';
%assignin('base','filename',filename)

% Clear all variables except filename
clear a* b* c* d* e* g* h* i* j* k* l* m* n* o* p* q* r* s* t* u* v* w* x* y* z*
clear A* B* C* D* E* G* H* I* J* K* L* M* N* O* P* Q* R* S* T* U* V* w* X* Y* Z*


% Parse filename from user input, and then 
% cd to directory where file will be saved to see if it already exists
find_filename_01 = strfind(filename,'/');
find_filename_02 = strfind(filename,'.log');
find_filename_03 = filename(find_filename_01(end)+1:find_filename_02-1);
pdf_file = [find_filename_03,'.pdf'];
cd_str = ['cd /F/Data/Analyzed/Level_I/',...
        filename(find_filename_01(end-1)+1:find_filename_01(end-1)+6)];
assignin('base','cd_str',cd_str)
eval(cd_str)
LS = ls;

if isempty(findstr(LS,pdf_file)) == 1 
else
    warndlg(['The file ',pdf_file,...
            ' already exists.  Please manually remove the file from the directory and re-execute the function.'],...
            'File Already Exists!')
break
end
        



% Open a log file to read
fid_in = fopen(filename,'r');

% Create a counter for the while loop
count = 0;

% Each iteration of the while loop reads a line from
% the file until the loop reaches the end of the file
while ~feof(fid_in)
    
inline = fgetl(fid_in);


    if ~isempty(inline)
    
        
        
        if inline(1) == '['
            find_colon = findstr(inline,':');
            find_brackets = findstr(inline,'[]');
            userinpt = inline(find_colon(1)+8:find_brackets-1);
            find_Trial = findstr(inline,'Trial :');
            find_curly_braces = findstr(inline,'{');
            sequence = inline(find_curly_braces+1:find_Trial-3);
            find_slash = findstr(inline,'/');
            date = [inline(find_slash(1)-1:find_slash(1)+7)];
            time = [inline(find_colon(1)-3:find_colon(1)+5)];

                if  inline(3) == ']' 
                %if isempty(inline(2)) == 1
                    trial_no = count+1;
                    %trial_no = str2num(inline(2));
                    log_cell{trial_no,1} = trial_no;
                    log_cell{trial_no,2} = userinpt;
                    log_cell{trial_no,3} = sequence;
                    log_cell{trial_no,14} = date;
                    log_cell{trial_no,15} = time;
                elseif inline(4) == ']'
                    %trial_no = str2num([inline(2) inline(3)]);
                    trial_no = count + 1;
                    log_cell{trial_no,1} = trial_no;
                    log_cell{trial_no,2} = userinpt;
                    log_cell{trial_no,3} = sequence;
                    log_cell{trial_no,14} = date;
                    log_cell{trial_no,15} = time;
                elseif inline(5) == ']'
                    %trial_no = str2num([inline(2) inline(3) inline(4)]);
                    trial_no = count+1;
                    log_cell{trial_no,1} = trial_no;
                    log_cell{trial_no,2} = userinpt;
                    log_cell{trial_no,3} = sequence;
                    log_cell{trial_no,14} = date;
                    log_cell{trial_no,15} = time;
                elseif inline(6) == ']'
                    trial_no = str2num([inline(2) inline(3) inline(4) inline(5)]);
                    trial_no = count+1;
                    log_cell{trial_no,1} = trial_no;
                    log_cell{trial_no,2} = userinpt;
                    log_cell{trial_no,3} = sequence;
                    log_cell{trial_no,14} = date;
                    log_cell{trial_no,15} = time;
                    
                %continue
                end
                  
                count = count+1;
                
               
        elseif isempty(strfind(inline,'program')) == 0
                
                 trial_no = count;
                 program = inline;
                 log_cell{trial_no,4} = program;
                 

        elseif isempty(strfind(inline,'delay')) == 0 & ...
                    isempty(strfind(inline,'pulse width')) == 0 & ...
                    isempty(strfind(inline,'voltage')) == 0
                 trial_no = count;
                 delay = inline(5:21);
                 pulse_width = inline(26:48);
                 voltage = inline(53:71);
                 log_cell{trial_no,5} = delay;
                 log_cell{trial_no,6} = pulse_width;
                 log_cell{trial_no,7} = voltage;
                 
        elseif isempty(strfind(inline,'num pulses')) == 0
                 trial_no = count;
                 num_pulses = inline;
                 log_cell{trial_no,8} = num_pulses;      
                
        elseif isempty(strfind(inline,'stimulator')) == 0
                 trial_no = count;
                 stimulator = inline;
                 log_cell{trial_no,9} = stimulator;
        elseif isempty(strfind(inline,'acquisition channels')) == 0 & ...
                   isempty(strfind(inline,'samples')) == 0 & ...
                   isempty(strfind(inline,'Offset')) == 0
                   trial_no = count;
                   acq_channels = inline(5:31);
                   samples = inline(34:46);
                   sample_rate = inline(51:64);
                   offset = inline(68:99);
                   log_cell{trial_no,10} = acq_channels;
                   log_cell{trial_no,11} = samples;
                   log_cell{trial_no,12} = sample_rate;
                   log_cell{trial_no,13} = offset;
        end
    end
end

fclose(fid_in)



% Label fields
log_cell{size(log_cell,1)+1,1} = 0;
log_cell(2:size(log_cell,1),:) = ...
     log_cell(1:end-1,:);
log_cell{1,1} = 'Trial No.';
log_cell{1,2} = 'User Inpt';
log_cell{1,3} = 'Sequence';
log_cell{1,4} = 'Program';
log_cell{1,5} = 'Delay';
log_cell{1,6} = 'Pulse W';
log_cell{1,7} = 'Voltage';
log_cell{1,8} = 'No. pulses';
log_cell{1,9} = 'Hardware';
log_cell{1,10} = 'Channels';
log_cell{1,11} = 'No. samples';
log_cell{1,12} = 'Samp rate';
log_cell{1,13} = 'Trig Offset';
log_cell{1,14} = 'Date';
log_cell{1,15} = 'Time';

% % % Replace empty cells with appropriate values if necessary
% If no program ...
for i = 2:size(log_cell,1)
    if isempty(log_cell{i,4}) == 1
       log_cell{i,4} = ' program : none';
    end
end

% If no delay ...
for i = 2:size(log_cell,1)
    if isempty(log_cell{i,5}) == 1
       log_cell{i,5} = 'delay : 0.0000';
    end
end

% If no pulse width ...
for i = 2:size(log_cell,1)
    if isempty(log_cell{i,6}) == 1
       log_cell{i,6} = 'pulse width :    0.0000';
    end
end


% If no voltage ...
for i = 2:size(log_cell,1)
    if isempty(log_cell{i,7}) == 1
       log_cell{i,7} = 'voltage : 0.0000';
    end
end

% If no num pulses ...
for i = 2:size(log_cell,1)
    if isempty(log_cell{i,8}) == 1
       log_cell{i,8} = '    num pulses : 0';
    end
end

% If no hardware ...
for i = 2:size(log_cell,1)
    if isempty(log_cell{i,9}) == 1
       log_cell{i,9} = 'none';
    end
end

% Ask if continuation of another recording
same_neuron_quest_box = 'Continuation';
same_neuron_answer_1 = questdlg(['Is this file a continuation of a log file from a previously recorded neuron?'],...
    same_neuron_quest_box,'No','Yes','No');
if strcmp(same_neuron_answer_1,'No') == 1
else
    same_neuron_answer_2_title = 'Continuation';
    same_neuron_answer_2_parameters = {'Which neuron is this a continuation of?'};
    same_neuron_answer_2 = inputdlg(same_neuron_answer_2_parameters,same_neuron_answer_2_title);
end



% Check channels present
channels = 0;
chan_1 = char(log_cell{2:end,10});
chan_2 = chan_1(:,23:26);
chan_3 = unique(chan_2,'rows');
for i = 1:size(chan_3,1)
    chan_str_str_1 = ['chan_str_',num2str(i),...
            ' = str2num(strrep(chan_3(',num2str(i),',:),''..'','':''));'];
    eval(chan_str_str_1)
    chan_str_str_2 = ['channels(1,size(channels,2)+1:size(channels,2)+size(chan_str_',num2str(i),',2)) = chan_str_',num2str(i),';'];
    eval(chan_str_str_2)
end
%assignin('base','chan_str_1',chan_str_1)
find_channels = find(channels);
channels = channels(1,find_channels);
channels = unique(channels)';
channels_str_1 = num2str(channels);
channels_str_2 = ['Channels ',channels_str_1(1),' through ',channels_str_1(end), ' were detected in log file.'];
channel_quest_box = 'Channel(s) Detected';
channel_answer_1 = questdlg([channels_str_2,' Were the sources of any of the channels changed during acquisition?'],...
    channel_quest_box,'No','Yes','No');
if channel_answer_1 == 'No'
    for i = 1:size(channels,1)
        Label_channels_parameters_str = ['Label_channels_parameters{1,',num2str(i),'} = [''Channel ',channels_str_1(i),'''];'];
        eval(Label_channels_parameters_str)
    end
    Label_channels_dlgbox_title = 'Label Channels';
    Label_channels_userinpt = inputdlg(Label_channels_parameters,Label_channels_dlgbox_title,'OK');
    %assignin('base','Label_channels_userinpt',Label_channels_userinpt)
end

for i = 1:size(Label_channels_userinpt,1)
    if isempty(Label_channels_userinpt{i})
        Label_channels_userinpt{i} = '<none>';
    end
end

% for i = 1:size(Label_channels_userinpt,1)
%     Label_channels_userinpt_str_str = ['Label_channels_userinpt_str = Label_channels_userinpt{',num2str(i),',1};'];
%     eval(Label_channels_userinpt_str_str)
%     assignin('base','Label_channels_userinpt_str',Label_channels_userinpt_str)
% %     if strcmp(Label_channels_userinpt_str,'DCN') == 1
% %         Label_channels_userinpt_str = 'DCN';
% %     else
% %     end
% end
% Separate conditions based on differences in trial parameters
Trial_params_char = [char(log_cell(2:end,2)),char(log_cell(2:end,3)),char(log_cell(2:end,4)),...
                     char(log_cell(2:end,5)),char(log_cell(2:end,6)),char(log_cell(2:end,7)),...
                     char(log_cell(2:end,8)),char(log_cell(2:end,9)),char(log_cell(2:end,10)),...
                     char(log_cell(2:end,11)),char(log_cell(2:end,12)),char(log_cell(2:end,13)),...
                     char(log_cell(2:end,14))];
Trial_params_cell =  cellstr(Trial_params_char);
%Trial_params_cell_uni = unique(Trial_params_cell,'rows')
%  for i = 2:size(log_cell,1)
%       str_comp_Trial_params_cell_str = ['str_comp_Trial_params_cell_Tr_',num2str(i-1),' = strcmp(Trial_params_cell('...
%             num2str(i-1),',1),Trial_params_cell(1:end,1));'];
%       eval(str_comp_Trial_params_cell_str)
%       Trial_params_mat_str = ['Trial_params_mat(:,',num2str(i-1),...
%               ') = str_comp_Trial_params_cell_Tr_',num2str(i-1),';'];
%       eval(Trial_params_mat_str)
%  end
 

[uni_el,last_uni_el,ind_uni_el] = unique(Trial_params_cell);
last_uni_el = sort(last_uni_el);
size_uni_el = size(uni_el,1);
for i = 1:size_uni_el
    % Groups the unique conditions together
    cond_str_str_1 = ['cond_',num2str(i),' = find(ind_uni_el == i);'];
    eval(cond_str_str_1)
end
    
% Place conditions in a matrix to be sorted
for i = 1:size_uni_el
cond_str_str_2 = ['cond_mat{1,',num2str(i),'} = cond_',num2str(i),';'];
eval(cond_str_str_2)
end
for i = 1:size(cond_mat,2)
    cond_mat_sort_str = ['cond_mat_sort(1,',num2str(i),') = cond_mat{1,',num2str(i),'}(1);'];
    eval(cond_mat_sort_str)
end
[cond_mat_sort,cond_mat_sort_i] = sort(cond_mat_sort); 
for i = 1:size(cond_mat_sort_i,2)
    cond_str_str_3 = ['cond_',num2str(i+size_uni_el),' = cond_',num2str(cond_mat_sort_i(i)),';'];
    eval(cond_str_str_3)
end
for i = 1:size_uni_el
    cond_str_str_4 = ['cond_',num2str(i),' = cond_',num2str(i+size_uni_el),';'];
    eval(cond_str_str_4)
end

for i = 1:size_uni_el   
    % Find the difference between successive indexes (or trials) of the
    % same condition
    cond_diff_str = ['cond_diff_',num2str(i),' = diff(cond_',num2str(i),');'];
    eval(cond_diff_str)
    %assignin('base','cond_diff_1',cond_diff_1)
    % Make sure the differences are unique (incomplete as of now)
    cond_diff_unique_str = ['cond_diff_uni_',num2str(i),...
            ' = unique(cond_diff_',num2str(i),');'];
    eval(cond_diff_unique_str)
    %assignin('base','cond_diff_uni_1',cond_diff_uni_1)
    % Create a string of the trial numbers for each unique condition
    if isempty(eval(['cond_diff_uni_',num2str(i)])) == 1
       cond_str_str = ['cond_str_',num2str(i),' = [num2str(cond_',num2str(i),')];'];
       eval(cond_str_str)
    else
        cond_str_str = ['cond_str_',num2str(i),' = [num2str(cond_',num2str(i),'(1)),',''':''',',num2str(cond_diff_uni_',num2str(i),'),',''':''',',num2str(cond_',num2str(i),'(end))];'];
        %assignin('base','cond_str_str',cond_str_str)
        eval(cond_str_str)
    end
end
%assignin('base','cond_str_1',cond_str_1)

% Dialog for unique conditions
if size_uni_el == 1
    Label_uni_cond_str_1 = [num2str(size_uni_el),' unique condition was found.  Label'...
        ' the condition in the next dialog box.'];
Label_uni_cond_str_2 = 'Unique Condition Detected.';
else
Label_uni_cond_str_1 = [num2str(size_uni_el),' unique conditions were found.  Label'...
        ' each unique condition in the next dialog box.'];
Label_uni_cond_str_2 = 'Unique Conditions Detected.';
end

Label_uni_cond_dlg_1 = questdlg(Label_uni_cond_str_1,Label_uni_cond_str_2,'OK','Cancel','OK');

% For each unique condition, create a string with 
% the appropriate userinput and sequence to be used in the label unique
% condition dialog.  Place these strings in a cell.
for i = 1:size(last_uni_el,1)
    Cond_Uni_str_str = ['cond_uni_',num2str(i),...
            '_str = [log_cell{last_uni_el(',num2str(i),...
            ')+1,2} log_cell{last_uni_el(',num2str(i),')+1,3}];'];
    eval(Cond_Uni_str_str)
end

for i = 1:size(last_uni_el,1)
    Cond_Uni_cell_str = ['cond_uni_cell{',num2str(i),...
            ',1} = cond_uni_',num2str(i),'_str;'];
    eval(Cond_Uni_cell_str)
end

% Create Label_unique_condition dialog box sized appropriately
% for the number of conditions found
if size(cond_uni_cell,1) <=10
    Label_uni_cond_str_3 = 'Label Unique Conditions.';
    Label_uni_cond_dlg_2 = inputdlg(cond_uni_cell,Label_uni_cond_str_3,1,cond_uni_cell);
else
    Label_uni_cond = figure;
    set(gcf,'Position',[200 200 1000 400]);
for i = 1:size(cond_uni_cell,1)
    if 10 >= i >= 1
        edit_box_str = ['edit_box_',num2str(i),' = uicontrol(gcf,''Style'',''edit'',''Position'',[10 405-40*',...
            num2str(i),' 400 20],''String'',''default'',''BackgroundColor'',[1 1 1],''String'',cond_uni_cell(i));'];
        eval(edit_box_str)
%         assignin('base','edit_box_str',edit_box_str)
%         assignin('base','edit_box_1',edit_box_1)
        Label_uni_cond_dlg_2_str = ['Label_uni_cond_dlg_2{',...
            num2str(i),',1} = char(get(edit_box_',num2str(i),',''String''));'];
        eval(Label_uni_cond_dlg_2_str)
        assignin_str = ['assignin(''base'',''edit_box_',num2str(i),''',edit_box_',num2str(i),')'];
        eval(assignin_str)
        assignin('base','Label_uni_cond_dlg_2',Label_uni_cond_dlg_2)
    else
        edit_box_str = ['edit_box_',num2str(i),' = uicontrol(gcf,''Style'',''edit'',''Position'',[450 40*',...
                num2str(i),'-75-80*(i-11) 400 20],''BackgroundColor'',[1 1 1],''String'',cond_uni_cell(i));'];
        eval(edit_box_str)
        Label_uni_cond_dlg_2_str = ['Label_uni_cond_dlg_2{',...
            num2str(i),',1} = char(get(edit_box_',num2str(i),',''String''));'];
        eval(Label_uni_cond_dlg_2_str)
        assignin_str = ['assignin(''base'',''edit_box_',num2str(i),''',edit_box_',num2str(i),')'];
        eval(assignin_str)
        assignin('base','Label_uni_cond_dlg_2',Label_uni_cond_dlg_2)
    end
end

assignin('base','Label_uni_cond',Label_uni_cond)
% Update Label_uni_cond_dlg_2 with any user changes
Label_uni_cond_button = uicontrol(gcf,'Style','pushbutton','Position',[900 350 50 30],...
    'String','OK','Callback',...
    ['for i = 1:size(Label_uni_cond_dlg_2,1),'...
     'edit_box_str_2 = [''edit_box_'',num2str(i),''_revised = get(edit_box_'',num2str(i),'',''''String'''');''];',...
     'eval(edit_box_str_2);',...
     'Label_uni_cond_dlg_2_str_2 = [''Label_uni_cond_dlg_2{'',num2str(i),'',1} = char(edit_box_'',num2str(i),''_revised);''];',...
     'eval(Label_uni_cond_dlg_2_str_2);',...
     'end; close(Label_uni_cond)']);
waitfor(Label_uni_cond)

end

% Find online user comments

fid_in = fopen(filename,'r');
FID_IN = fscanf(fid_in,'%c');
find_userinpt_brackets_1 = findstr(FID_IN,'[]');
FID_IN(find_userinpt_brackets_1) = [];
find_userinpt_brackets_2 = findstr(FID_IN,'] {');
FID_IN(find_userinpt_brackets_2) = [];
find_com_begin = findstr(FID_IN,'/*');
find_com_end = findstr(FID_IN,'*/');
find_trace_open_brackets = findstr(FID_IN,'[');
find_trace_close_brackets = findstr(FID_IN,']');

for i = 1:size(find_com_begin,2)/2
    online_user_com_begin_str = ['online_user_com_begin_',num2str(i),' = find_com_begin(2*',num2str(i),'-1);'];
    eval(online_user_com_begin_str)
    online_user_com_trace_linspace_str = ['online_user_com_trace_linspace_',num2str(i),' = linspace(online_user_com_begin_',...
            num2str(i),',online_user_com_begin_',num2str(i),',size(find_trace_open_brackets,2));'];
    eval(online_user_com_trace_linspace_str)
    online_user_com_trace_min_diff_str = ['[y,online_user_com_trace_min_diff_',...
        num2str(i),'] = min(abs(find_trace_open_brackets - online_user_com_trace_linspace_',num2str(i),'));'];
    eval(online_user_com_trace_min_diff_str)
    online_user_com_trace_str = ['online_user_com_trace_',...
            num2str(i),' = FID_IN(find_trace_open_brackets(online_user_com_trace_min_diff_',num2str(i),')+1:find_trace_close_brackets(online_user_com_trace_min_diff_',num2str(i),')-1);'];
    eval(online_user_com_trace_str)
    online_user_com_str_str = ['online_user_com_str_',num2str(i),' = FID_IN(find_com_begin(2*',num2str(i),'-1)+129:find_com_end(',num2str(i),'*2)-128);'];
    eval(online_user_com_str_str)
end


% % Write .com file
% Parse filename from user input
find_filename_01 = strfind(filename,'/');
find_filename_02 = strfind(filename,'.log');
find_filename_03 = filename(find_filename_01(end)+1:find_filename_02-1);
com_file = [find_filename_03,'_com'];
ps_file = [find_filename_03,'.ps'];
pdf_file = [find_filename_03,'.pdf'];

% cd to directory where file will be saved
cd_str = ['cd /F/Data/Analyzed/Level_I/',...
        filename(find_filename_01(end-1)+1:find_filename_01(end-1)+6)];
assignin('base','cd_str',cd_str)
eval(cd_str)

% Begin writing file
fid = fopen(com_file,'w');
text_asterisks = '******************************';
text_filename = find_filename_03;
text_channels = 'Channel(s) Recorded:';

fprintf(fid,'%s \s\s',text_asterisks);

if exist('same_neuron_answer_2') == 1
    fprintf(fid,'%s \s\s',text_filename);
    fprintf(fid,'%s \n\n',['(same neuron as ',char(same_neuron_answer_2),')']);
else
    fprintf(fid,'%s \n',text_filename);
end
fprintf(fid,'%s \t',text_channels);

for i = 1:size(Label_channels_userinpt,1)
    if i == 1
        fprintf(fid,'%s \n\t\t\t',[num2str(channels(i)),' - ',char(Label_channels_userinpt(i))]);
    elseif i == size(Label_channels_userinpt,1)
        fprintf(fid,'%s \n\n\n',[num2str(channels(i)),' - ',char(Label_channels_userinpt(i))]);
    else
        fprintf(fid,'%s \n\t\t\t',[num2str(channels(i)),' - ',char(Label_channels_userinpt(i))]);
    end
end

% Check for all channels being same
All_channels = char(log_cell{2:end,10});
All_channels_uni = unique(All_channels,'rows')
if size(All_channels_uni,1) == 1
    fprintf(fid,'%s \n\n','All Channels:')
end

for i = 1:size(Label_uni_cond_dlg_2,1)
    size_cond_str_str = ['size_cond_str = size(cond_str_',num2str(i),',2);'];
    eval(size_cond_str_str)
    if size_cond_str == 1
        trials_userinpt_str = ['fprintf(fid,''%s \n\n'',[cond_str_',num2str(i),',''               ''',',char(Label_uni_cond_dlg_2(i))]);'];
        eval(trials_userinpt_str)
    elseif size_cond_str == 2
        trials_userinpt_str = ['fprintf(fid,''%s \n\n'',[cond_str_',num2str(i),',''              ''',',char(Label_uni_cond_dlg_2(i))]);'];
        eval(trials_userinpt_str)
    elseif size_cond_str == 6
        trials_userinpt_str = ['fprintf(fid,''%s \n\n'',[cond_str_',num2str(i),',''          ''',',char(Label_uni_cond_dlg_2(i))]);'];
        eval(trials_userinpt_str)
    elseif size_cond_str == 7
        trials_userinpt_str = ['fprintf(fid,''%s \n\n'',[cond_str_',num2str(i),',''         ''',',char(Label_uni_cond_dlg_2(i))]);'];
        eval(trials_userinpt_str)
    elseif size_cond_str == 8
        trials_userinpt_str = ['fprintf(fid,''%s \n\n'',[cond_str_',num2str(i),',''        ''',',char(Label_uni_cond_dlg_2(i))]);'];
        eval(trials_userinpt_str)
    else
        trials_userinpt_str = ['fprintf(fid,''%s \n\n'',[cond_str_',num2str(i),',''       ''',',char(Label_uni_cond_dlg_2(i))]);'];
        eval(trials_userinpt_str)
    end
end

fprintf(fid,'%s \n\n',text_asterisks);

if exist('online_user_com_str_1') == 1
    fprintf(fid,'%s \n\n','**** Online User Comments ****');
                     
    for i = 1:size(find_com_begin,2)/2
        user_com_trace_str = ['fprintf(fid,''%s \n\n'',[''Trace#: '',online_user_com_trace_',num2str(i),']);'];
        eval(user_com_trace_str)
        user_com_str_str = ['fprintf(fid,''%s \n\n'',['''''''',online_user_com_str_',num2str(i),','''''''']);'];
        eval(user_com_str_str)
    end

    fprintf(fid,'%s \n\n\n',text_asterisks);
end

% Ask if user wants to save, convert to ps and then to a pdf and display using acroread
save_quest_box = 'Save';
save_answer_1 = questdlg(['Would you like to save this file?'],save_quest_box,'Yes','No','Yes');
if strcmp(save_answer_1,'Yes') == 1
    
    % Add offline comments
    add_offline_quest_box = 'Offline Comments';
    add_offline_com_answer = questdlg(['Would you like to add offline comments?'],add_offline_quest_box,'Yes','No','Yes');
    if strcmp(add_offline_com_answer,'Yes') == 1
        offline_com_parameters = {'Add offline comments below:'};
        offline_com_title = 'Offline comments';
        offline_com_userinpt = inputdlg(offline_com_parameters,offline_com_title,6);
        fprintf(fid,'%s \n\n','**** Offline User Comments ***');
        for i = 1:size(char(offline_com_userinpt),1)
            fprintf(fid,'%s \n\n',char(offline_com_userinpt{1}(i,:)));
        end
        fprintf(fid,'%s \n\n\n',text_asterisks);
    else
    end
     
    % Add spike discrimination parameters
    spike_discrim_quest_box = 'Spike Discrimination Parameters';
    spike_discrim_answer_1 = questdlg(['Would you like to save spike discrimination parameters?'],...
    spike_discrim_quest_box,'Yes','No','Yes');
    if strcmp(spike_discrim_answer_1,'Yes') == 1
        Label_channels_spike_discrim_params_dlgbox_title = 'Label Spike Discrimination Parameters';
        Label_channels_spike_discrim_params_userinpt = inputdlg(Label_channels_userinpt,Label_channels_spike_discrim_params_dlgbox_title,'OK');
        for i = 1:size(Label_channels_userinpt,1)
            if isempty(Label_channels_spike_discrim_params_userinpt{i,1})
                Label_channels_and_spike_discrim_params_str = ['Label_channels_and_spike_discrim_params{',num2str(i),',1} = [num2str(channels(',num2str(i),'))  '' - '' char(Label_channels_userinpt(',num2str(i),'))];'];
                eval(Label_channels_and_spike_discrim_params_str)
            else
                Label_channels_and_spike_discrim_params_str = ['Label_channels_and_spike_discrim_params{',num2str(i),',1} = [num2str(channels(',num2str(i),')) '' - '' char(Label_channels_userinpt(',num2str(i),')),'' - comments: '',char(Label_channels_spike_discrim_params_userinpt(',num2str(i),'))];'];
                eval(Label_channels_and_spike_discrim_params_str)
            end
        end
            
%                 spike_discrim_title = 'Spike Discrimination Parameters';
%                 spike_discrim_parameters = {'Min Amp','Max Amp','Min Width','Max Width'};
%                 spike_discrim_userinpt = inputdlg(spike_discrim_parameters,spike_discrim_title);
%                 sed_str = ['!sed ''s/DCN/DCN  ',...
%                     ['(discrimination parameters: ',char(spike_discrim_userinpt{1}),',',char(spike_discrim_userinpt{2}),',',...
%                     char(spike_discrim_userinpt{3}),',',char(spike_discrim_userinpt{4}),')'],'/'' ',com_file,' > com_file2'];
%                 eval(sed_str)
%         mv_com_file_str = ['!mv ',com_file,' com_file2'];
%         eval(mv_com_file_str)
%     end   
%     else
        mv_com_file_str = ['!mv ',com_file,' com_file2'];
        eval(mv_com_file_str)
    else
        mv_com_file_str = ['!mv ',com_file,' com_file2'];
        eval(mv_com_file_str)
    end
    
    fclose(fid)
    ps_str = ['!a2ps -B -R --columns=1 com_file2 -o ',ps_file];
    eval(ps_str)
    ps_str = ['!ps2pdfwr ',ps_file,' ',pdf_file];
    eval(ps_str)
    acroread_str = ['!acroread ',pdf_file,'&'];
    eval(acroread_str)
    ps_str_2 = ['!rm ',ps_file];
    eval(ps_str_2)
    com_str = ['!rm ',com_file];
    eval(com_str)
    com_str_2 = ['!rm com_file2'];
    eval(com_str_2)
    
else
    fclose(fid)
    ps_str = ['!a2ps -B -R --columns=1 ',com_file,' -o ',ps_file];
    eval(ps_str)
    ps_str = ['!ps2pdfwr ',ps_file,' ',pdf_file];
    eval(ps_str)
    ps_str_2 = ['!rm ',ps_file];
    eval(ps_str_2)
    com_str = ['!rm ',com_file];
    eval(com_str)
    acroread_str = ['!acroread ',pdf_file,'&'];
    eval(acroread_str)
    pause(5)
    delete_pdf_str = ['!rm ',pdf_file];
    eval(delete_pdf_str)
end


% cd back to /F/Data/Raw 
cd_str_2 = ['cd /F/Data/Raw'];
eval(cd_str_2)





