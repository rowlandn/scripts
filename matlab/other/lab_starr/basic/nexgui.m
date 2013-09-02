function nexgui()
%Loads the main gui and sets things up

clear all;
global h_gui % handle made invisible so it is not corrupted by command line users, 
             % but need to get at it sometimes, so make it global
global NEXROOT; % RLM: location of THIS Nex source file

%set(0, 'FixedWidthFont', 'Tahoma');
fprintf('Building GUI, this may take a few seconds...');
me = mfilename('fullpath'); % RLM: to fetch code directory
[NEXROOT,name,ext,versn] = fileparts([me '.m']); % RLM: save for later!

if ~isempty(findobj('Tag','NexGUI'))
    h_gui = openfig('nexgui.fig','reuse');
else
    h_gui = openfig('nexgui.fig','reuse');
    movegui(h_gui,'northeast'); % RLM: TOP-RIGHT corner
    set(h_gui,'HandleVisibility','on','Units','characters');
    % initialize here
    set(h_gui,'HandleVisibility','Callback'); % RLM: protects it from commandline
end
fprintf(' ready.\n');

format short g; %number format of command window
