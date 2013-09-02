function nexguifunctions(action)
%
%

switch action
    case 'Chk_MergeNexfiles'
        if get(gco,'Value')
            set(findobj('Tag','Edt_NumNexMerged'),'Enable','on');
        else
            set(findobj('Tag','Edt_NumNexMerged'),'Enable','off');
        end
        
    case 'Edt_NumNexMerged'
        val = str2num(get(gco,'String'));
        if val ~= round(val)
            set(gco,'String',num2str(val));
        end
        if isempty(val) | (val<2)
            set(gco,'String',num2str(2));
        end
        
    case 'Chk_AutomateFileloads'
        
        
    case 'Chk_ISICalculations'
        if get(gco,'Value')
            set(findobj('Tag','Chk_ISIUseThreshold'),'Enable','on','Value',0);
            set(findobj('Tag','Edt_ISIMinimumThreshold'),'Enable','on');
            set(findobj('Tag','Chk_ISIBatching'),'Enable','on','Value',0);
        else
            set(findobj('Tag','Chk_ISIUseThreshold'),'Enable','off','Value',0);
            set(findobj('Tag','Edt_ISIMinimumThreshold'),'Enable','off');
            set(findobj('Tag','Chk_ISIBatching'),'Enable','off','Value',0);
        end
        
    case 'Chk_ISIUseThreshold'
        if get(gco,'Value')
            set(findobj('Tag','Edt_ISIMinimumThreshold'),'Enable','on');
        else
            set(findobj('Tag','Edt_ISIMinimumThreshold'),'Enable','off');
        end
        
    case 'Edt_ISIMinimumThreshold'
        val = str2num(get(gco,'String'));
        if isempty(val) | (val<0.1)
            set(gco,'String',num2str(1));
        end
        
    case 'Chk_ISIBatching'
        if get(gco,'Value')
            set(findobj('Tag','Chk_ISICalculations'),'Enable','off');
            set(findobj('Tag','Chk_MergeNexfiles'),'Enable','off');
            set(findobj('Tag','Chk_AutomateFileloads'),'Enable','off');
            set(findobj('Tag','Chk_MergeContinuous'),'Enable','off');
            set(findobj('Tag','Chk_AutomateVariableLoads'),'Enable','off');
            set(findobj('Tag','Pop_AutomateFrom'),'Enable','off');
            set(findobj('Tag','Chk_MergeCogentLog'),'Enable','off');
            set(findobj('Tag','Chk_MergeCogentLogKeystrokes'),'Enable','off');
            set(findobj('Tag','Chk_WriteNex'),'Enable','off');
            set(findobj('Tag','Chk_RenameComment'),'Enable','off');
            set(findobj('Tag','Chk_RenameVariables'),'Enable','off');
        else
            nexguifunctions('default')
%            set(gco,'Enable','on');
%             set(findobj('Tag','Chk_ISICalculations'),'Enable','on');
%             set(findobj('Tag','Chk_MergeNexfiles'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_AutomateFileloads'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_MergeContinuous'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_AutomateVariableLoads'),'Enable','on','Value',0);
%             set(findobj('Tag','Pop_AutomateFrom'),'Enable','on');
%             set(findobj('Tag','Chk_MergeCogentLog'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_MergeCogentLogKeystrokes'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_WriteNex'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_RenameComment'),'Enable','on','Value',0);
%             set(findobj('Tag','Chk_RenameVariables'),'Enable','on','Value',0);
        end
        
    case 'Chk_MergeContinuous'
        if get(gco,'Value')
            set(findobj('Tag','Chk_AutomateVariableLoads'),'Enable','on','Value',0);
        else
            set(findobj('Tag','Chk_AutomateVariableLoads'),'Enable','off','Value',0);
            set(findobj('Tag','Pop_AutomateFrom'),'Enable','off');
        end
        
    case 'Chk_AutomateVariableLoads'
        if get(gco,'Value')
            set(findobj('Tag','Pop_AutomateFrom'),'Enable','on');
        else
            set(findobj('Tag','Pop_AutomateFrom'),'Enable','off');
        end
        
    case 'Pop_AutomateFrom'
        
    case 'Chk_MergeCogentLog'
        if get(gco,'Value')
            set(findobj('Tag','Chk_MergeCogentLogKeystrokes'),'Enable','on');
        else
            set(findobj('Tag','Chk_MergeCogentLogKeystrokes'),'Enable','off');
        end
        
    case 'Chk_MergeCogentLogKeystrokes'
        
    case 'Chk_WriteNex'
        if get(gco,'Value')
            set(findobj('Tag','Chk_RenameComment'),'Enable','on');
            set(findobj('Tag','Chk_RenameVariables'),'Enable','on');
        else
            set(findobj('Tag','Chk_RenameComment'),'Enable','off');
            set(findobj('Tag','Chk_RenameVariables'),'Enable','off');
        end
        
    case 'Chk_RenameComment'
        
    case 'Chk_RenameVariables'
        
    case 'default' % DEFAULT
            set(findobj('Tag','Chk_MergeNexfiles'),'Enable','on','Value',0);
            set(findobj('Tag','Edt_NumNexMerged'),'Enable','off');
            set(findobj('Tag','Chk_AutomateFileloads'),'Enable','on','Value',1);
            set(findobj('Tag','Chk_ISICalculations'),'Enable','on','Value',0);
            set(findobj('Tag','Chk_ISIUseThreshold'),'Enable','off','Value',0);
            set(findobj('Tag','Edt_ISIMinimumThreshold'),'Enable','off');
            set(findobj('Tag','Chk_ISIBatching'),'Enable','off');            
            set(findobj('Tag','Chk_MergeContinuous'),'Enable','on','Value',1);
            set(findobj('Tag','Chk_AutomateVariableLoads'),'Enable','off','Value',0);
            set(findobj('Tag','Pop_AutomateFrom'),'Enable','off');
            set(findobj('Tag','Chk_MergeCogentLog'),'Enable','on','Value',0);
            set(findobj('Tag','Chk_MergeCogentLogKeystrokes'),'Enable','off','Value',0);
            set(findobj('Tag','Chk_WriteNex'),'Enable','on','Value',1);
            set(findobj('Tag','Chk_RenameComment'),'Enable','on','Value',0);
            set(findobj('Tag','Chk_RenameVariables'),'Enable','on','Value',0);
            
    otherwise
            % just in case
end
        
        
        
        
        