function separate_color_traces(Selected_Trace_no)

Temp_filename_raw = evalin('base','Temp_filename_raw');
max_Selected_Trace_no = max(Temp_filename_raw(:,Selected_Trace_no))
min_Selected_Trace_no = min(Temp_filename_raw(:,Selected_Trace_no))

%[sort_Traces,sort_Traces_i] = sort(max_Selected_Trace_no)






 colors = ['b','r','g','c','m','k','y'];
 
%figure
% for i = 1:size(Selected_Trace_no,2)
%     if i == 1
%         hold on
%         plot(Temp_filename_raw(:,Selected_Trace_no(i)),colors(i))
%     elseif i >= 2 & i<=7
% plot(Temp_filename_raw(:,Selected_Trace_no(i)),colors(i))




for i = 1:size(Selected_Trace_no,2)
    if i == 1
        hold on
        plot(Temp_filename_raw(:,Selected_Trace_no(i)),colors(i))
    elseif i >= 2 & i<=7
        min_last_Selected_Trace_no(i) = min(Temp_filename_raw(:,Selected_Trace_no(i-1)))
        max_Selected_Trace_no = max(Temp_filename_raw(:,Selected_Trace_no(i)))
        plot(Temp_filename_raw(:,Selected_Trace_no(i))-abs((max_Selected_Trace_no-sum(min_last_Selected_Trace_no(1:i)))),colors(i))
        %plot(Temp_filename_raw(:,Selected_Trace_no(i)),colors(i))
    end
end

% for i = 1:size(Selected_Trace_no,2)
%     legend_parameters{1,i} = ['Trace ',num2str((Selected_Trace_no(i)))];
% end
% 
% legend_dialogbox_title = 'Enter Legend';
% legend_userinpt = inputdlg(legend_parameters,legend_dialogbox_title);
% legend(legend_userinpt)

