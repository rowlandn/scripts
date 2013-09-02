function varargout = SpikeRaster(varargin)
% SpikeRaster()
% This program reads in NEX sorted spikes and generates a raster plot
% Created by RST, 2003-06-25
%
%      H = SPIKERASTER returns the handle to a new SPIKERASTER or the handle to
%      the existing singleton*.
%
%      SPIKERASTER('Property','Value',...) creates a new SPIKERASTER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SpikeRaster_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SPIKERASTER('CALLBACK') and SPIKERASTER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SPIKERASTER.M with the given input
%      arguments.

% Last Modified by GUIDE v2.5 26-Jun-2003 17:12:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpikeRaster_OpeningFcn, ...
                   'gui_OutputFcn',  @SpikeRaster_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before SpikeRaster is made visible.
function SpikeRaster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

handles.tickcolor = 'k';		% ticks are black='k' by default

[nexname, RTPATH] = uigetfile('*.nex', 'Select spike time file (NEX)');

if (nexname ~= 0)
    cd(RTPATH);
    [nvar, varname, types] = nex_info(nexname);
    if nvar > 1
        error([num2str(nvar) ' spikes in ' nexname '!!!  I don''t know how to process > 1 spike yet']);
    end
    [spk.n, spk.t] = nex_ts(nexname,varname);
else
    error(['I can''t find the NEX file:  ' nexname ' in ' RTPATH]);
end

% convert spk.t to spike times in msec.s
handles.spk_t = round(1000*spk.t);        % put spk time in msec

% Set default values
handles.width = 1000;
handles.length_of_raster = handles.spk_t(end);
handles.start = 0;
handles.tickwidth = 0.05;
handles.tickheight = 0.8;

% Update fname displayed on gui
set(handles.current_file,'String', nexname);

% Choose default command line output for SpikeRaster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpikeRaster wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpikeRaster_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function width_Callback(hObject, eventdata, handles)
handles.width = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function length_of_raster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_of_raster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function length_of_raster_Callback(hObject, eventdata, handles)
handles.length_of_raster = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes on button press in file_open.
function file_open_Callback(hObject, eventdata, handles)
[nexname, RTPATH] = uigetfile('*.nex', 'Select spike time file (NEX)');

if (nexname ~= 0)
    cd(RTPATH);
    [nvar, varname, types] = nex_info(nexname);
    if nvar > 1
        error([num2str(nvar) ' spikes in ' nexname '!!!  I don''t know how to process > 1 spike yet']);
    end
    [spk.n, spk.t] = nex_ts(nexname,varname);
else
    error(['I can''t find the NEX file:  ' nexname ' in ' RTPATH]);
end

% convert spk.t to spike times in msec.s
handles.spk_t = round(1000*spk.t);        % put spk time in msec

% Update fname displayed on gui
set(handles.current_file,'String', nexname);
guidata(hObject, handles);


% --- Executes on button press in plot_raster.
function plot_raster_Callback(hObject, eventdata, handles)

num_spaces = 1-handles.tickheight;

% Calculation of start of displayed data strem
begin = threshold(handles.spk_t,handles.start,'-+');
if begin == -1,begin = 1;end;

% Calculation of end of displayed data stream
finish = threshold(handles.spk_t,handles.start+handles.length_of_raster,'-+');
if finish == -1,finish = length(handles.spk_t)+1;end;

% Display of data stream
D = handles.spk_t(begin:(finish-1));% Selection of data segment to be displayed
D = D-handles.start;						% Adjustment so that the data segment starts with 0
X = rem(D,handles.width);							% Calculation of X-axis of plot
Y = floor(D/handles.width)*(num_spaces+handles.tickheight);	% Calculation of Y-axis of plot, adjusting for line spacing

figure
plotbar_color(X,Y,handles.tickheight,handles.tickwidth,handles.tickcolor);		% raster plot modified version of Scott's function
axis tight
axis ij				% put axis in image mode

% fix Y axis to reflect time in seconds
v = str2num(get(gca,'YTickLabel')) * handles.width/1000;	
set(gca,'YTickLabel',num2str(v));

set(gca,'TickDir','out');
title(get(handles.current_file,'String'),'Interpreter','none');

% --- Executes during object creation, after setting all properties.
function start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function start_Callback(hObject, eventdata, handles)
handles.start = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tickwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tickwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function tickwidth_Callback(hObject, eventdata, handles)
handles.tickwidth = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tickheight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tickheight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tickheight_Callback(hObject, eventdata, handles)
handles.tickheight = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over start.
function start_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


