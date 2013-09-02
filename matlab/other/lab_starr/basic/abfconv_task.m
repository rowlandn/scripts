%  abfconv_task()%% This function imports an .abf file collected by the Axon%   Assumes a sampling rate in file of 20 kHz per channel% %   Channel 1 is converted into a two binary files for spike sorting%       1) "fname.wav" = a file which can be played  -> MSD sorting (AlphaOmega)%       2) "fname.raw" = a file which can be imported into Offline Sorter (Plexon)% RLM:  3) "fname.ddt" = a ddt file which is native to Plexon, for Offline Sorter% RLM:  4) "fname-accel.ddt" = an OPTIONAL ddt file containing accelerometer for accel deflection "sorting"%%   Channels 2-n (assummed to contain EMG data) are saved in "fname.mat" for subsequent merging w/ spike times%       EMG data are processed as follows:%           1) rectified%           2) low pass filtered%           3) sub-sampled -> 1 kHz (1 sample/msec)%%   The spike channel is also saved in "fname.mat" as a 1 kHz "envelope" representation%       (max & min values across 2 msec are stored in successive msec samples)%% RLM: created by Rob Turner, last edited 2003.10.03 by Rory Michaelis% Revised to not filter channel 5 - contains event dataclear all% RLM: ##################################################################rename_variables = 1 % RLM: TOGGLE prompt for manually naming emg channelswrite_accel_ddt = 1  % RLM: TOGGLE writing of LAST analog channel as accelerometer DDTaccel_channel_append = '_Accel'; % RLM: default name appendage for Accelerometer channel (in *.MAT)accel_ddt_append = '_a'; % RLM: name appendage for accelerometer DDT (no extension!!!)% RLM: ##################################################################fn_dir = strrep(which('abfconv2'),'abfconv2.m','');   % Store function location for later usedisp('**CHOOSE AN .ABF FILE AND MAKE SURE TO HIGHLIGHT ALL DESIRED CHANNELS BEFORE IMPORTING');i=1;while i    ABFML % choose abf file from CD or proper directory    % Check to see if spike channel imported (AD0)    if and( size(who('*AD0')) > 0,  size(who('*Time')) > 0)        break;    end    if size(who('*AD0')) > 0        ButtonName = questdlg('Spike chan (AD0) not selected! Try again?');        switch ButtonName,        case 'No',            break        case 'Cancel',            return        end % switch    end    if size(who('*Time')) > 0        ButtonName = questdlg('Time not imported! Try again?');        switch ButtonName,        case 'No',            break        case 'Cancel',            return        end % switch    endend% Compute sampling rate of abf file (it is not always the same)if size(who('*Time')) > 0	tmp = cell2struct(who('*Time'),'name',1);	time = eval(tmp.name); % let time = [TIME array of the .abf file)    fs = round(length(time)/time(end));         % Sampling rate    RATE1 = fs/4000;  % Output EMG is first sub-sampled at 1/RATE1    RATE2 = 4000/1000;  % Output EMG is sub-sampled later at 1/RATE2	RATE_PLT = fs/1000;else    ButtonName = questdlg('Time was not imported! Assuming 50 microsec. sampling interval.');    fs = 20000;         % Sampling rate    RATE = 20;  % Output EMG is sub-sampled later at 1/RATEend% Load filter informationfilt_file = 'abfconv_filters';load([ fn_dir filt_file ]);% If Spike channel importedif size(who('*AD0')) > 0	tmp = cell2struct(who('*AD0'),'name',1);	ad0 = eval(tmp.name); % let ad0 = [ADO array of the .abf file)	outputname = strrep(tmp.name,'_1_AD0','');	    % Filter out HF (~9 kHz) noise common in ABF files    flt_unit = filtfilt(LowPass_5kHz_20kHz.tf.num, 1, ad0 );        % Rescale data to an appropriate range	if abs(min(flt_unit)) > max(flt_unit)         rng = abs(min(flt_unit)); % set range to largest amplitude	else         rng = max(flt_unit);	end	rng = rng*1.1;      % add a little to range to ensure we don't clip output data		%flt_unit component has the spike data	scaledad0 = (32767/rng)*flt_unit;%     fid = fopen ([outputname '.raw'],'wb');	% 	disp(['Writing spikes channel to:  ' outputname '.raw']);% 	fwrite (fid,scaledad0,'int16');% 	fclose (fid);% %     disp(['Writing spikes channel to:  ' outputname '.wav']);%     wavwrite((1/rng)*flt_unit,fs,16,[outputname '.wav']);% RLM: ##################################################################	ddtdata = (32767/(max(abs(flt_unit))))*flt_unit; % +/-2047 value is from DDT format 12bit res.        writeddt(ddtdata,fs,outputname,16,'Creating ddt file.');% RLM: ##################################################################        clear(tmp.name);    % Finally, remove Spk channel from memoryend% Now process AUX chans if imported	% get AUX channel names and process    tmp = cell2struct(who('*AD*'),'name',2);    [n_aux x] = size(tmp);    ndata = length(eval(tmp(1).name));    	% Make decimated plot to help evt channel selection	plt_inds = 1:round(RATE_PLT):ndata;	h = figure;	for i =1:n_aux		subplot(n_aux,1,i)		eval([ 'plot(time(plt_inds),' tmp(i).name '(plt_inds), ''-k'');']);		title(['Chan' num2str(i)]);	end		% Select Event channel	evt_ch = 0;	while evt_ch < 1 | evt_ch > n_aux		evt_ch = input('Enter channel # for Task Events: ');	end	emg_ch = [];	emg_ch = input('Enter channel #s (if any) for good EMG (format: [1 2 ...]): ');	accel_ch = input('Enter channel # for ACCEL: ');	if exist('h','var')		if ishandle(h)			close( h )		end	end	    % For first down sampling from ~20kHz -> 4kHz    ndeci1 = round(ndata/RATE1);       % Length of decimated data    deci_list1 = zeros(1,ndeci1);    for i=1:ndeci1                    % Make deci_list the hard way to make sure we subsample         deci_list1(i) = round(i*RATE1);%   correctly for unusual sampling rates    end    if deci_list1(end) > length(time)    % Avoid overrun because of rounding        deci_list1(end) = length(time);    end        % For 2nd down-sampling from 4kHz -> 1 kHz    ndata1 = length(deci_list1);    deci_list2 = 1:RATE2:ndata1;        aux_chan = zeros(n_aux,length(deci_list2));        fprintf('Processing  %d  EMG/Event/Accel channels', n_aux);	for i = 1:n_aux        ch = eval(tmp(i).name);        if i ~=accel_ch & i~=evt_ch            ch = detrend(ch);					% Remove offset & linear trend            ch = abs(ch);						% Compute absolute value		end		fprintf('.');		if i~=evt_ch			ch = filtfilt(LowPass_2kHz_20kHz.tf.num, 1, ch ); % Fist filter out > 2 kHz 		end		deci1 = ch(deci_list1);					% Sub-sample down to 4 kHz		fprintf('.');		if i~=evt_ch			deci1 = filtfilt(LowPass_100Hz_4kHz.tf.num, 1, deci1 ); % Filter out > 100 Hz 		end		aux_chan(i,:) = deci1(deci_list2);		% Sub-sample down to 1 kHz		fprintf('.');	end	fprintf('\n');    % Include time vector if imported	tmp = cell2struct(who('*Time'),'name',1);	time = eval(tmp.name);	time = squeeze(time(deci_list1));   % Now sub-sample	time = squeeze(time(deci_list2));   % Now sub-sample	% Detect events	events = ParseIntraopEvents(time,aux_chan(evt_ch,:));    if ~isempty(events) & ~isempty(accel_ch) & accel_ch>0 & accel_ch<=n_aux        events = DetectAccels(time,aux_chan(accel_ch,:),events);    end		    outputname = [outputname '.mat'];	disp(['Writing EMG/ACCEL  channels to:  ' outputname]);	save( outputname, 'n_aux', 'aux_chan', 'time', 'events');