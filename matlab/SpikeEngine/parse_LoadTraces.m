function parse_LoadTraces(LoadTraces_userinpt)


% clear; clc
% LoadTraces_userinpt{1,1} = '1';
% LoadTraces_userinpt{2,1} = '1:4:197'
% if isempty(findstr(LoadTraces_userinpt{2},'-')) == 0
%     warndlg('Please use the colon separator in place of dashes.','Warning!')
% else
% Check for spaces.  The number of arrays in the input equals the number
% of spaces in the input plus one.
LoadTraces_userinpt_spaces = findstr(char(LoadTraces_userinpt(2)),' ');
size_LoadTraces_userinpt = LoadTraces_userinpt_spaces + 1;

% If there are no spaces, accept user input
if isempty(LoadTraces_userinpt_spaces) == 1
    LoadTraces_num_str = ...
        ['LoadTraces_num = str2num(char(LoadTraces_userinpt{2,1}));'];
    eval(LoadTraces_num_str)
    LoadTraces_str = num2str(LoadTraces_num);
    assignin('base','LoadTraces_str',LoadTraces_str)
else

% If there spaces, identify individual arrays
for i = 1:size(LoadTraces_userinpt_spaces,2) + 1
    if i == 1 & diff([i LoadTraces_userinpt_spaces(i)]) == 1
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),' = LoadTraces_userinpt{2}(',...
            num2str(i),')'];
        eval(LoadTrace_str_str)
    elseif i == 1 & diff([i LoadTraces_userinpt_spaces(i)]) ~= 1
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),' = LoadTraces_userinpt{2}(',...
            num2str(i),':LoadTraces_userinpt_spaces(',...
            num2str(i),')-1)'];
        eval(LoadTrace_str_str)
    elseif i == size(LoadTraces_userinpt_spaces,2) & ...
            diff([LoadTraces_userinpt_spaces(i-1)+1 LoadTraces_userinpt_spaces(i)]) == 1
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),...
        ' = LoadTraces_userinpt{2}(LoadTraces_userinpt_spaces(',num2str(i),')-1)'];
        eval(LoadTrace_str_str)
    elseif i == size(LoadTraces_userinpt_spaces,2) & ...
            diff([LoadTraces_userinpt_spaces(i-1)+1 LoadTraces_userinpt_spaces(i)]) ~= 1
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),...
        ' = LoadTraces_userinpt{2}(LoadTraces_userinpt_spaces(',num2str(i-1),')+1:LoadTraces_userinpt_spaces(',num2str(i),')-1)'];
        eval(LoadTrace_str_str)
    elseif i ~= size(LoadTraces_userinpt_spaces,2)+1 & diff([i+1 LoadTraces_userinpt_spaces(i)]) == 1
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),...
                ' = LoadTraces_userinpt{2}(LoadTraces_userinpt_spaces(',num2str(i-1),')+1)'];
        eval(LoadTrace_str_str)
    elseif i ~= size(LoadTraces_userinpt_spaces,2)+1 & diff([i+1 LoadTraces_userinpt_spaces(i)]) ~= 1
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),' = LoadTraces_userinpt{2}(',...
        '(LoadTraces_userinpt_spaces(',num2str(i-1),...
        ')+1):(LoadTraces_userinpt_spaces(',num2str(i),')-1))'];
        eval(LoadTrace_str_str)
    elseif i == size(LoadTraces_userinpt_spaces,2)+1 
        LoadTrace_str_str = ['LoadTrace_str_',num2str(i),...
                ' = LoadTraces_userinpt{2}(LoadTraces_userinpt_spaces(',num2str(i-1),...
                ')+1:size(LoadTraces_userinpt{2},2))'];
        eval(LoadTrace_str_str)
    end
end


%Place individual array strings into a cell.  Check for:
% - valid array
% - positive integers
 
%Places individual array string into a cell
for i = 1:size(LoadTraces_userinpt_spaces,2) + 1
    LoadTrace_cell_str_str = ['LoadTrace_cell_str{i,1} = LoadTrace_str_',num2str(i)];
    eval(LoadTrace_cell_str_str)
    % Checks for positive numbers
    if isempty(findstr(char(LoadTrace_cell_str(i,:)),'-')) == 0 
        warndlg(['''',char(LoadTrace_cell_str(i)),''' is an invalid trace number. Please enter positive integers for trace numbers.'],'Warning!')
    break    
    end
    % Checks for integers
    if isempty(findstr(char(LoadTrace_cell_str(i,:)),'.')) == 0 
        warndlg(['''',char(LoadTrace_cell_str(i)),''' is an invalid trace number. Please enter positive integers for trace numbers.'],'Warning!')
    break
    end
    %Checks for valid array and then coverts string to numbers
    if isempty(findstr(char(LoadTrace_cell_str(i,:)),':')) == 0 
        find_colon = findstr(char(LoadTrace_cell_str(i,1)),':')
        if size(find_colon,2) == 1
            LoadTrace_cell_str_num_beg = ['LoadTrace_cell_num_beg = str2num(LoadTrace_cell_str{',num2str(i),',1}(1:find_colon-1))'];
            eval(LoadTrace_cell_str_num_beg)
            LoadTrace_cell_str_num_end = ['LoadTrace_cell_num_end = str2num(LoadTrace_cell_str{',num2str(i),',1}(find_colon+1:end))'];
            eval(LoadTrace_cell_str_num_end) 
                if LoadTrace_cell_num_beg > LoadTrace_cell_num_end
                    warndlg(['''',char(LoadTrace_cell_str(i)),''' is an invalid array.'],'Warning!')
                break
                else
                    LoadTrace_str_num = ['LoadTrace_num_',num2str(i),' = str2num(LoadTrace_str_',...
                    num2str(i),')'];
                    eval(LoadTrace_str_num)
                end
        elseif size(find_colon,2) == 2
            LoadTrace_cell_str_num_beg = ['LoadTrace_cell_num_beg = str2num(LoadTrace_cell_str{',num2str(i),',1}(1:find_colon(1)-1))'];
            eval(LoadTrace_cell_str_num_beg)
            LoadTrace_cell_str_num_end = ['LoadTrace_cell_num_end = str2num(LoadTrace_cell_str{',num2str(i),',1}(find_colon(2)+1:end))'];
            eval(LoadTrace_cell_str_num_end) 
            if LoadTrace_cell_num_beg > LoadTrace_cell_num_end
                    warndlg(['''',char(LoadTrace_cell_str(i)),''' is an invalid array.'],'Warning!')
                    break
            else
                LoadTrace_str_num = ['LoadTrace_num_',num2str(i),' = str2num(LoadTrace_str_',...
                    num2str(i),')'];
                eval(LoadTrace_str_num)
            end
        end
    else
    LoadTrace_str_num = ['LoadTrace_num_',num2str(i),' = str2num(LoadTrace_str_',...
            num2str(i),')'];
    eval(LoadTrace_str_num)
    end
    
end

% Place numbers into single array
for i = 1:size(LoadTraces_userinpt_spaces,2) + 1
    if i == 1
        LoadTrace_num_str = ['LoadTraces_num(1,1:size(LoadTrace_num_',num2str(i),',2)) = LoadTrace_num_',num2str(i)];
    else
        LoadTrace_num_str = ['LoadTraces_num(1,size(LoadTraces_num,2)+1:size(LoadTraces_num,2)+size(LoadTrace_num_',num2str(i),',2)) = LoadTrace_num_',num2str(i)];
    end
    eval(LoadTrace_num_str)
end

% Eliminate duplicate traces and sort final array
LoadTraces_num_sort = sort(LoadTraces_num)
diff_LoadTraces_num_sort = diff(LoadTraces_num_sort)
find_zero_diff_LoadTraces_num_sort = find(diff_LoadTraces_num_sort == 0)
%LoadTraces_num_sort(1,find_zero_diff_LoadTraces_num_sort) =[]
LoadTraces_num = LoadTraces_num_sort;
LoadTraces_str = num2str(LoadTraces_num);
end

assignin('base','LoadTraces_num',LoadTraces_num)
assignin('base','LoadTraces_str',LoadTraces_str)

% Create Trace_nos_cell for SpikeEngine
no_traces = size(LoadTraces_num,2)
for i = 1:size(LoadTraces_num,2)
    if size(num2str(LoadTraces_num(i)),2) == 1
        Trace_nos_cell_str = ['Trace_nos_cell(',num2str(i),',:)= [''  Raw_Trace                Trace_0'',num2str(LoadTraces_num(i)),''                    Channel_'',num2str(char(LoadTraces_userinpt{3}))''];'];
        eval(Trace_nos_cell_str)  
    elseif size(num2str(LoadTraces_num(i)),2) == 2
        Trace_nos_cell_str = ['Trace_nos_cell(',num2str(i),',:)= [''  Raw_Trace                Trace_'',num2str(LoadTraces_num(i)),''                    Channel_'',num2str(char(LoadTraces_userinpt{3}))''];'];
        eval(Trace_nos_cell_str)  
    end
end
assignin('base','Trace_nos_cell',Trace_nos_cell)