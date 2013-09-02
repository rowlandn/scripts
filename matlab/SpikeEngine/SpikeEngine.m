%[right up width height]

% Clear all variables, close all figures and clear the command window
clear all; close all; clc

% Analysis_mode_question_box = 'Data Analysis Mode';
% Analysis_mode_answer = questdlg('Select Data Analysis Mode',Analysis_mode_question_box,'Manual','Automatic','Manual');
% if Analysis_mode_answer == 'Manual'
%     Parse_log_question_box = 'Parse Data Log';
%     Parse_log_answer = questdlg('Would you like to parse the data log file before beginning analysis?',Parse_log_question_box,'Yes','No','Yes');
%     if strcmp(Parse_log_answer,'Yes')
%         [pcdx_log_filename,pcdx_log_pathname, pcdx_log_filterindex] = uigetfile('*.log', 'Select a .log file');
%         parse_pcdx_log([pcdx_log_pathname,pcdx_log_filename])
%     else
%     end
% else
% end

% Display splashscreen
load splash_sound
spl = figure;
set(gcf,'Color',[1 1 1],'Position',[380 278 500 300],'ToolBar','none','Menu','none');
set(gca,'Box','off','Color',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1])
text(0,.4,'Loading Burst-Triggered Averaging...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075,'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Spike-Triggered Averaging...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Trace Averaging...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Autocorrelation...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Crosscorrelation...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Burst Detection (L-S)...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Response Detection (GLR)...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading ISI Histogram...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Spectrogram...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Coherence...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Partial Coherence...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Direct Transfer Function...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Make Movie...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Parse PCDX...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Parse NeuroSage...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Record Macro...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Filter Traces...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Invert Traces...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Resample Signal...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Subtract Artifact...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
cla
text(0,.4,'Loading Subtract mean...');
t1 = text(0,0,'Nathan C. Rowland, Dieter Jaeger, Alfonso Delgado, Jeremy Edgerton, Tom Sangrey,');
t2 = text(0,-.075, 'Cengiz Gunay, Su Li');
set(t1,'FontSize',8)
set(t2,'FontSize',8)
pause(.15)
%cla
%pause(1)
sound(splash_sound,22050)
close


% Create figure (and axis) with the following specifications
fig = figure
screen_size = get(0,'ScreenSize');
set(gcf,'Position',[15 40 1250 700],'MenuBar','none');
subplot(10,10,[1:7,11:17,21:27,31:37,41:47,51:57,61:67,71:77,81:87])
axis_2
zoom_SE

subplot(10,10,[1:7,11:17,21:27,31:37,41:47,51:57,61:67,71:77,81:87])

% Just for aesthetics
Traces_text = uicontrol(gcf,'Style','text','Position',[900 230 300 420],'String',['Trace Type      Trace Number        Channel Number'],'FontWeight','bold','FontSize',10,'BackgroundColor',[.8 .8 .8]);
Prev_button = uicontrol(gcf,'Style','pushbutton','Position',[900 210 30 25],'String','<<');
Forward_button = uicontrol(gcf,'Style','pushbutton','Position',[930 210 30 25],'String','>>');
Add_button = uicontrol(gcf,'Style','pushbutton','Position',[900 185 30 25],'String','+');
Remove_button = uicontrol(gcf,'Style','pushbutton','Position',[930 185 30 25],'String','-');
Select_all_button = uicontrol(gcf,'Style','pushbutton','Position',[1140 210 60 25],'String','Select All');
Separate_traces_button = uicontrol(gcf,'Style','checkbox','Position',[1020 210 120 25],'String','Separate Traces');
Trace_no_listbox = uicontrol(gcf,'Style','listbox','Max',2,'Position', [900 235 300 400],'BackgroundColor',[1 1 1]);
Graphs_text = uicontrol(gcf,'Style','text','Position',[1000 175 90 15],'String','Graphs:','FontWeight','bold','FontSize',10,'BackgroundColor',[.8 .8 .8]);
Graph_str_listbox = uicontrol(gcf,'Style','listbox','Position', [900 75 300 100],'BackgroundColor',[1 1 1]);
     


% Create title
if isunix == 1
    if screen_size(4) == 1024
        Title = uicontrol(gcf,'Style','text','Position',[390 770 200 30], 'String',...
        'Spike Engine','FontSize',20,'BackgroundColor',[.8 .8 .8]);
    elseif screen_size(4) == 800
        Title = uicontrol(gcf,'Style','text','Position',[390 650 200 30], 'String',...
        'Spike Engine','FontSize',20,'BackgroundColor',[.8 .8 .8]);
    end
else
    if screen_size(4) == 800
        Title = uicontrol(gcf,'Style','text','Position',[390 660 200 30], 'String',...
        'Spike Engine','FontSize',20,'BackgroundColor',[.8 .8 .8]);
    end 
end

% Menu - Advanced Plot
SpikeEngine_Menu_Advanced_Plot

% Menu - Analysis
SpikeEngine_Menu_Analysis

% Menu - Clear
SpikeEngine_Menu_Clear

% Menu - Tools
SpikeEngine_Menu_Tools

% Menu - Transform
SpikeEngine_Menu_Transform



% Function - LoadTraces 
SpikeEngine_Function_LoadTraces

% Function - FindSpikes 
SpikeEngine_Function_FindSpikes
     
% Function - PlayTraces 
SpikeEngine_Function_PlayTraces

% Function - GenerateGraphs 
SpikeEngine_Function_GenerateGraphs