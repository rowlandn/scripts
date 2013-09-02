for i = 1:size(LoadTraces_num,2)
    if size(num2str(LoadTraces_num(i)),2) == 1
        Spike_Trace_nos_cell_str = ['Spike_Trace_nos_cell(',num2str(i),',:)= [''  Spike_Trace                Trace_0'',num2str(LoadTraces_num(i)),''                    Channel_'',num2str(char(LoadTraces_userinpt{3}))''];'];
        eval(Spike_Trace_nos_cell_str)  
    elseif size(num2str(LoadTraces_num(i)),2) == 2
        Spike_Trace_nos_cell_str = ['Spike_Trace_nos_cell(',num2str(i),',:)= [''  Spike_Trace                Trace_'',num2str(LoadTraces_num(i)),''                    Channel_'',num2str(char(LoadTraces_userinpt{3}))''];'];
        eval(Spike_Trace_nos_cell_str)  
    end
end

assignin('base','Spike_Trace_nos_cell',Spike_Trace_nos_cell)