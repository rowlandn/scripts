%function abf_spk()
% abf_spk()
% This function imports an .abf file 
%   Corrects sampling rate according to time array returned by abfml
% 	Read matching nex file (spike times)
%	Plots raw data
%	Plots windowed spike data around spike times in nex file
%

clear all
fn_dir = strrep(which('abf_spk'),'abf_spk.m','');   % Store function location for later use


disp('**CHOOSE AN .ABF FILE AND MAKE SURE TO HIGHLIGHT ALL DESIRED CHANNELS BEFORE IMPORTING');
i=1;
while i
    abfml % choose abf file from CD or proper directory

    % Check to see if spike channel imported (AD0)
    if and( size(who('*AD0')) > 0,  size(who('*Time')) > 0)
        break;
    end
    if size(who('*AD0')) > 0
        ButtonName = questdlg('Spike chan (AD0) not selected! Try again?');
        switch ButtonName,
        case 'No',
            break
        case 'Cancel',
            return
        end % switch
    end
    if size(who('*Time')) > 0
        ButtonName = questdlg('Time not imported! Try again?');
        switch ButtonName,
        case 'No',
            break
        case 'Cancel',
            return
        end % switch
    end
end

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


% Compute sampling rate of abf file (it is not always the same)
if size(who('*Time')) > 0
	tmp = cell2struct(who('*Time'),'name',1);
	time = eval(tmp.name); % let time = [TIME array of the .abf file)
    fs = round(length(time)/time(end));         % Sampling rate
	smp_interv = 1/fs;
    RATE1 = fs/4000;  % Output EMG is first sub-sampled at 1/RATE1
    RATE2 = 4000/1000;  % Output EMG is sub-sampled later at 1/RATE2
else
    ButtonName = questdlg('Time was not imported! Assuming 50 microsec. sampling interval.');
    fs = 20000;         % Sampling rate
    RATE = 20;  % Output EMG is sub-sampled later at 1/RATE
end

% Load filter information
filt_file = 'abfconv_filters';
load([ fn_dir filt_file ]);

% If Spike channel imported
if size(who('*AD0')) > 0
	tmp = cell2struct(who('*AD0'),'name',1);
	ad0 = eval(tmp.name); % let ad0 = [ADO array of the .abf file)
	outputname = strrep(tmp.name,'_1_AD0','');
	
    % Filter out HF (~9 kHz) noise common in ABF files
    flt_unit = filtfilt(LowPass_5kHz_20kHz.tf.num, 1, ad0 );
end

% Plot data for this many spikes
min_spk = 100;
max_spk = 200;
% Define time window in msec pre- post-spike trigger
spk_time_win = -0.6:smp_interv*1000:0.8;
% Compute indices for spike window
spk_ind_win = round(spk_time_win(1)*fs/1000) : round(spk_time_win(end)*fs/1000);

figure
hold on
spk.ind = zeros(max_spk-min_spk,1);
for i = min_spk:max_spk
	spklo = spk.t(i) - smp_interv/2;
	spkhi = spk.t(i) + smp_interv/2;
	spk.ind(i-min_spk+1) = find( time>spklo & time<spkhi );
	
	plot(spk_time_win,flt_unit(spk_ind_win + spk.ind(i-min_spk+1)),'-k');
end
xlabel('msec');
ylabel('volts');

% Make plot of continuous data during same period as sorted spikes
contin_ind_win = (spk_ind_win(1)+spk.ind(1)):(spk_ind_win(end)+spk.ind(end));
figure
plot(time(contin_ind_win),flt_unit(contin_ind_win),'-k');
xlabel('msec');
ylabel('volts');
