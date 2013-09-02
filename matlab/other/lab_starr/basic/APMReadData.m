function apmdata = APMReadData( filename, split, substract_noise )
% apmdata = APMReadUserData( filename, split )
%  Reads data from an APM data file
%  Input arguments:
%   filename    - input data file
%   split       - a flag that when set to 1 (true), data is split in trials every
%                 time a start of trial code is encountered in the data
%                 file. Typically, this flag is set to 1 when recording APM data
%                 in gated mode with timers reset.
%   substract_noise - when set to 1 (true), the digital noise synchronized with 
%                 USB packet transmission is removed using an adaptive filtering algorithm
%
%  Output:
%   apmdata - a structure containing all the data
%          in the file. Type "apmdata" at matlab command prompt to 
%          list all the memebers. When using gated recording mode, the
%          function returns a vector of structures, each trial being 
%          stored in a different element of the vector.
%
%

% Message codes definitions, see FileMessageCodes.h file in the C API
MSG_NACC_WAVEFORM           = 00;
MSG_ACC_WAVEFORM	        = 01;
MSG_NOTRIAL_NACC_WAVEFORM	= 02;
MSG_NOTRIAL_ACC_WAVEFORM	= 03;
MSG_TIMESTAMP				= 04;
MSG_START_OF_TRIAL			= 05;
MSG_END_OF_TRIAL			= 06;
MSG_EXT_EVENT				= 07;
MSG_TRIGGER_LEVEL			= 08;
MSG_TRIGGER_TIME			= 09;
MSG_TRIGGER_SLOPE			= 10;
MSG_WINDOW_TIME				= 11;
MSG_WINDOW_LOW				= 12;
MSG_WINDOW_HIGH				= 13;
MSG_TIME_CALIBRATION		= 14;
MSG_VOLTAGE_CALIBRATION		= 15;
MSG_SAMPLING_FREQUENCY		= 16;
MSG_REWARD					= 17;
MSG_LFP_TIME_CALIBRATION	= 18;
MSG_DIGINPUT_CHANGED	    = 19;
MSG_STROBED_DATA        	= 20;
MSG_UART_DATA               = 21;
MSG_SPIKE_LEN       	    = 22;   % Code added by FSS
MSG_AVG_SAMPLES         	= 23;   % Code added by FSS
MSG_DRIVE_DEPTH            	= 24;
MSG_CHANNEL_INFO           	= 25;   % Type (if missing, then APM_CHAN is assumed)
MSG_CHANNEL_NAME           	= 26;   % Channel name, a string that's padded to a multiple of 4 characters, to fit into the message structure
MSG_LINE_NOISE_FILTER      	= 27;
MSG_HIGH_PASS_FREQ         	= 28;
MSG_LOW_PASS_FREQ         	= 29;


MSG_TEMPLATE				= 32;
MSG_CONTINUOUS				= 48;
MSG_LFP     				= 56;
MSG_SLOW_DATA     			= 64;
MSG_NOTRIAL_SLOW_DATA     	= 66;
MSG_AUX_DATA        	    = 72;
MSG_TCPIP_USER0				= hex2dec('10000');
MSG_TCPIP_USER1				= hex2dec('10001');     % Filename of the behavioral data for the current trial
MSG_TCPIP_USER2				= hex2dec('10002');
MSG_TCPIP_USER3				= hex2dec('10003');
MSG_TCPIP_USER4				= hex2dec('10004');
MSG_TCPIP_USER5				= hex2dec('10005');
MSG_TCPIP_USER6				= hex2dec('10006');
MSG_TCPIP_USER7				= hex2dec('10007');

% Channel types
APM_CHAN        = 5;
GL4KSP_CHAN     = 6;
GL4KSC_CHAN     = 7;
GL4KAUX_CHAN    = 8;

if (nargin < 1)
    error('The function requires at least one input argument ...\n   Please type "help APMRead" for a description of the function');
end;

if isempty(filename)
    [filename, pathname] = uigetfile('*.apm', 'APM Data File (*.apm)');
    filename = strcat(pathname, filename);
end;

try
    split;
catch
    split = 0;
end;

try
    substract_noise;
catch
    substract_noise = 0;        % Switch for substracting or not the interference with the digital transmissions
end;

% The unit is coded on bits 8 to 11 when using template matching 
%  or on bit 0 when using the acceptance window

% Initialization
nchan   = 128;              % Default number of channels
nunits  = 4;                % Default number of units per channel
trials  = [];
if (split)
    ntr = zeros(1,nchan);   % trial number; must keep track of it on each channel individually, since data may be interleaved between channels.
else
    ntr = ones(1,nchan);
end;
nts     = zeros(1,nchan);   % number of timestamps for each channel
last_ts = zeros(1,nchan);   % last timestamp for each channel
last_cts = zeros(1,nchan);   % last timestamp for continuous data on each channel
last_chts = zeros(1,nchan);   % last timestamp for chart(slow) data on each channel
nwf     = zeros(1,nchan);   % number of waveforms for each channel
nevt    = zeros(1,nchan);   % number of external events for each channel
nevt1   = zeros(1,nchan);   % number of reward/validate events for each channel
trig_level  = zeros(1,nchan);  % trigger level, in ADC levels
trig_slope  = zeros(1,nchan);  % trigger slope
trig_time   = zeros(1,nchan);  % trigger time, in ADC samples
win_time    = zeros(1,nchan);  % lower window level, in ADC levels
win_low     = zeros(1,nchan);  % upper window level, in ADC levels
win_high    = zeros(1,nchan);  % window position in time, in ADC samples
spike_len   = zeros(1,nchan);  % spike length, for internal use in Guideline4000
pkt_len     = zeros(1,nchan);  % packet length for continuous data
LNF         = zeros(1,nchan);  % line noise filter on/off
high_pass_freq = zeros(1,nchan);  % high-pass filter frequency, in Hz
low_pass_freq  = zeros(1,nchan);  % low-pass filter frequency, in Hz
tcal      = repmat(1e6/48000.0,1,nchan);  % default time calibration, in microseconds per ADC sample
lfptcal   = repmat(1e6/500.0,1,nchan);    % default time calibration for LFP data
charttcal = repmat(1e6/100.0,1,nchan);    % default time calibration for chart data
vcal    = repmat(1.41/32767.0,1,nchan);   % default voltage calibration
%chan_type    = repmat(APM_CHAN,1,nchan);  % default channel type (5=APM, 6=GL4KSP, 7=GL4KSC, 8=GL4KAUX)
chan_type    = repmat(GL4KSP_CHAN,1,nchan);  % default channel type (5=APM, 6=GL4KSP, 7=GL4KSC, 8=GL4KAUX)
chan_type(11) = GL4KAUX_CHAN;   % Test
chan_type(12) = GL4KAUX_CHAN;   % Test
sampfreq    = repmat(48000,1,nchan);  % default sampling frequency;  sampling frequency is not sent any more; use time calibration instead
current_filename = [];      % A user message sent through TCPIP, telling the filename of the current trial descriptor and data
nrlen = 512;                % This must be larger than the maximum packet length
min_cont_pkts = 10; % Minimum number of packets to calculate a noise replica
fixed_noise_freq = [4000 6000 8000];    % The fixed noise frequencies that may present in the signal, regardless of the other channels' sampling frequencies

fdat=fopen(filename,'rb');
if (fdat~=-1)
    disp(sprintf('Processing %s ...',filename));
    while (~feof(fdat))
        m_code=fread(fdat,1,'uint32');               % Message type
        if (~feof(fdat))
            m_channel=fread(fdat,1,'uint32');        % Message channel & unit
            m_unit=bitshift(m_channel,-16);          % Upper word of the channel may contain the unit for certain messages
            m_channel=bitand(m_channel,hex2dec('0000FFFF'));  % Lower word of the channel contains the channel
            m_length=fread(fdat,1,'uint32');         % Message length
            %disp([m_code m_channel m_length]);
            %disp(sprintf('Channel %d, message code %d, length %d',m_channel,m_code, m_length));
%             Deal here with all messages ...
            switch (m_code)
                case MSG_CHANNEL_INFO
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %disp(sprintf('Channel %d, type %d, firmware version %.2f !',m_channel,m_data(1),double(m_data(2))/100.0));
                    chan_type(m_channel)=m_data(1);
                case MSG_TIME_CALIBRATION
                    %disp('Time Calibration record ...');
                    tcal(m_channel)=fread(fdat,m_length,'float32');
                    sampfreq(m_channel)=round(1e6/tcal(m_channel));
                    disp(sprintf('Channel %d, time calibration %f ...',m_channel,tcal(m_channel)));
                case MSG_SAMPLING_FREQUENCY
                    %disp('Sampling frequency changed ...');
                    sampfreq(m_channel)=fread(fdat,m_length,'uint32');   % Read message data, a 32-bit integer
                    tcal(m_channel)=1e6/sampfreq(m_channel);
                    disp(sprintf('Channel %d, sampling frequency %d ...',m_channel,sampfreq(m_channel)));
                case MSG_START_OF_TRIAL
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    last_cts(m_channel)=0;
                    if (split)
                        ntr(m_channel)=ntr(m_channel)+1;
                    end;
                    if (chan_type(m_channel)==GL4KSP_CHAN)
                        trials(ntr(m_channel)).channels(m_channel).noise_replica = zeros(nrlen,1);      % Noise replica on each channel
                        trials(ntr(m_channel)).channels(m_channel).continuous_samples  = 0;  % Number of continuous ADC samples, for each channel, used for preallocating the data vectors (much faster processing)
                        trials(ntr(m_channel)).channels(m_channel).continuous_packets  = 0;  % Number of continuous ADC packets, for substracting the digital noise during recording
                    end;
                case MSG_CONTINUOUS
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %disp(sprintf('Channel %d, continuous data packet, length %d, timestamp %d',m_channel,m_length,m_data(1)));
                    if ((ntr(m_channel)>0) && (m_data(1)>last_cts(m_channel)) && (last_cts(m_channel)>0))
                        disp(sprintf(' %d samples dropped.',m_data(1)-last_cts(m_channel)));
                        % Pad the data vector with zeros, in the same number as the missing samples
                         trials(ntr(m_channel)).channels(m_channel).continuous_samples =  trials(ntr(m_channel)).channels(m_channel).continuous_samples + m_data(1)-last_cts(m_channel);
                    end;
                    pkt_len(m_channel) = max(pkt_len(m_channel), m_length-1);
                    last_cts(m_channel) = m_data(1)+m_length-1;   % This is in fact pointing to the expected timestamp of the next packet
                    if (ntr(m_channel)>0)
                        try
                            trials(ntr(m_channel)).channels(m_channel);
                        catch
                            % For compatibility with older file formats
                            trials(ntr(m_channel)).channels(m_channel).noise_replica = zeros(nrlen,1);      % Noise replica on each channel
                            trials(ntr(m_channel)).channels(m_channel).continuous_samples  = 0;  % Number of continuous ADC samples, for each channel, used for preallocating the data vectors (much faster processing)
                            trials(ntr(m_channel)).channels(m_channel).continuous_packets  = 0;  % Number of continuous ADC packets, for substracting the digital noise during recording
                        end;
                       
                        trials(ntr(m_channel)).channels(m_channel).noise_replica(1:(m_length-1)) = trials(ntr(m_channel)).channels(m_channel).noise_replica(1:(m_length-1)) + double(m_data(2:end));
                        trials(ntr(m_channel)).channels(m_channel).continuous_samples = trials(ntr(m_channel)).channels(m_channel).continuous_samples + m_length-1;
                        trials(ntr(m_channel)).channels(m_channel).continuous_packets = trials(ntr(m_channel)).channels(m_channel).continuous_packets + 1;
                    end;
                case MSG_AUX_DATA
                    %disp(sprintf('Auxiliary I/O data on channel %d',m_channel));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %m_data
                    %dec2hex(m_data)
                    %disp(sprintf('Aux data on channel %d, length %d, active channels 0x%s, timestamp %d ',m_channel,m_length,dec2hex(m_unit),m_data(1)));
                    pkt_len(m_channel) = max(pkt_len(m_channel), m_length-1);
                otherwise
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
            end;
        end;
    end;
    if (split==0)
        for i=1:nchan
            if (chan_type(i)==GL4KSP_CHAN)
                for j=1:ntr(i)
                    if (length(trials(j).channels)>=i)
                        %disp(sprintf('Channel %d, trial %d of %d',i,j,ntr(i)));
                        if (isfield(trials(j).channels(i),'continuous_packets'))
                            if (trials(j).channels(i).continuous_packets>min_cont_pkts)
                                trials(j).channels(i).noise_replica = trials(j).channels(i).noise_replica/trials(j).channels(i).continuous_packets;
                            end;
                            if (trials(j).channels(i).continuous_samples>0)
                                disp(sprintf('Channel %d, trial %d, continuous samples %d',i,j,trials(j).channels(i).continuous_samples));
                                trials(j).channels(i).continuous = zeros(trials(j).channels(i).continuous_samples,1);
                                trials(j).channels(i).continuous_samples = 0;   % Reset continuous_samples, as it will be re-used in the main loop
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    frewind(fdat);
    if (split)
        ntr = zeros(1,nchan);   % trial number; must keep track of it on each channel individually, since data may be interleaved between channels.
    else
        ntr = ones(1,nchan);
    end;
    last_cts = zeros(1,nchan);   % reset last timestamp for continuous data on each channel
    while (~feof(fdat))
        m_code=fread(fdat,1,'uint32');               % Message type
        if (~feof(fdat))
            m_channel=fread(fdat,1,'uint32');        % Message channel & unit
            m_unit=bitshift(m_channel,-16);          % Upper word of the channel may contain the unit for certain messages
            m_channel=bitand(m_channel,hex2dec('0000FFFF'));  % Lower word of the channel contains the channel
            m_length=fread(fdat,1,'uint32');         % Message length
            %disp([m_code m_channel m_length]);
            %disp(sprintf('Channel %d, message code %d, length %d',m_channel,m_code, m_length));
%             Deal here with all messages ...
            switch (m_code)
                case MSG_CHANNEL_INFO
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, type %d, firmware version %.2f !',m_channel,m_data(1),double(m_data(2))/100.0));
                    chan_type(m_channel)=m_data(1);
                case MSG_START_OF_TRIAL
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, Start-of-Trial at %d !',m_channel,m_data));
                    nts(m_channel)=0;
                    nwf(m_channel)=0;
                    last_ts(m_channel)=0;
                    last_cts(m_channel)=0;
                    nevt(m_channel)=0;      % Reset the event count
                    nevt1(m_channel)=0;      % Reset the event count
                    if (split)
                        ntr(m_channel)=ntr(m_channel)+1;
                    else
                        if (chan_type(m_channel)==GL4KAUX_CHAN)
                            try
                                trials(ntr(m_channel)).aux(m_channel).start_trial=[trials(ntr(m_channel)).aux(m_channel).start_trial m_data];
                            catch
                                trials(ntr(m_channel)).aux(m_channel).start_trial=m_data;
                            end;
                        else
                            try
                                trials(ntr(m_channel)).channels(m_channel).start_trial=[trials(ntr(m_channel)).channels(m_channel).start_trial m_data];
                            catch
                                trials(ntr(m_channel)).channels(m_channel).start_trial=m_data;
                            end;
                        end;
                    end;
                    if (chan_type(m_channel)~=GL4KAUX_CHAN)
                        % Set current trial parameters that are available so far ...
                        if (~isempty(current_filename))
                            trials(ntr(m_channel)).filename=current_filename;   % The current filename message is sent BEFORE the start-of-trial message
                        end;
                        trials(ntr(m_channel)).channels(m_channel).trig_level=trig_level(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).trig_time=trig_time(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).win_time=win_time(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).win_low=win_low(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).win_high=win_high(m_channel);
                        %trials(ntr(m_channel)).channels(m_channel).sampling_frequency=sampfreq(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).time_calibration=tcal(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).lfp_time_calibration=lfptcal(m_channel);
                        trials(ntr(m_channel)).channels(m_channel).voltage_calibration=vcal(m_channel);
                        for i=1:nunits
                            trials(ntr(m_channel)).channels(m_channel).template(i).unit=i;
                            try
                                trials(ntr(m_channel)).channels(m_channel).template(i).waveform=templates(m_channel).unit(i).waveform;
                            catch
                                trials(ntr(m_channel)).channels(m_channel).template(i).waveform=[];
                            end;
                        end;
                    end;
                case MSG_END_OF_TRIAL
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, End-of-Trial at %d !',m_channel,m_data));
                    if (~split)
                        if (chan_type(m_channel)==GL4KAUX_CHAN)
                            try
                                trials(ntr(m_channel)).aux(m_channel).end_trial=[trials(ntr(m_channel)).aux(m_channel).end_trial m_data];
                            catch
                                trials(ntr(m_channel)).aux(m_channel).end_trial=m_data;
                            end;
                        else
                            try
                                trials(ntr(m_channel)).channels(m_channel).end_trial=[trials(ntr(m_channel)).channels(m_channel).end_trial m_data];
                            catch
                                trials(ntr(m_channel)).channels(m_channel).end_trial=m_data;
                            end;
                        end;
                    end;
                case MSG_EXT_EVENT
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    disp(sprintf('External Event Input ON at %d',m_data));
                    if (ntr(m_channel)>0)
                        nevt(m_channel)=nevt(m_channel)+1;
                        %trials(ntr(m_channel)).channels(m_channel).events(nevt(m_channel))=m_data; % For APM-02 store event timestamp in the spike channel structure
                        trials(ntr(m_channel)).aux(m_channel).events(nevt(m_channel))=m_data; % For Guideline, store event timestamp in the auxiliary channel structure
                    end;
                case MSG_TRIGGER_LEVEL
                    trig_level(m_channel)=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, trigger level %d ...',m_channel,trig_level(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).trig_level=trig_level(m_channel);
                    end;
                case MSG_TRIGGER_SLOPE
                    trig_slope(m_channel)=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, trigger slope %d ...',m_channel,trig_slope(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).trig_slope=trig_slope(m_channel);
                    end;
                case MSG_TRIGGER_TIME
                    trig_time(m_channel)=fread(fdat,m_length,'uint32');     % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, trigger position %d ...',m_channel,trig_time(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).trig_time=trig_time(m_channel);
                    end;
                case MSG_WINDOW_TIME
                    win_time(m_channel)=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, window position in time %d ...',m_channel,win_time(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).win_time=win_time(m_channel);
                    end;
                case MSG_WINDOW_LOW
                    win_low(m_channel)=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, window low level %d ...',m_channel,win_low(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).win_low=win_low(m_channel);
                    end;
                case MSG_WINDOW_HIGH
                    win_high(m_channel)=fread(fdat,m_length,'int32');   % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, window high level %d ...',m_channel,win_high(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).win_high=win_high(m_channel);
                    end;
                case MSG_TIME_CALIBRATION
                    %disp('Time Calibration record ...');
                    tcal(m_channel)=fread(fdat,m_length,'float32');
                    disp(sprintf('Channel %d, time calibration %f ...',m_channel,tcal(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).time_calibration=tcal(m_channel);
                    end;
                case MSG_VOLTAGE_CALIBRATION
                    %disp('Voltage Calibration record ...');
                    vcal(m_channel)=fread(fdat,m_length,'float32');
                    disp(sprintf('Channel %d, voltage calibration %f ...',m_channel,vcal(m_channel)));
                    if (chan_type(m_channel)==GL4KAUX_CHAN)
                        if (ntr(m_channel)>0)
                            trials(ntr(m_channel)).aux(m_channel).voltage_calibration = vcal(m_channel);
                        end;
                    else
                        if (ntr(m_channel)>0)
                            trials(ntr(m_channel)).channels(m_channel).voltage_calibration=vcal(m_channel);
                        end;
                    end;
                case MSG_SAMPLING_FREQUENCY         % These messages are not sent any more. Use time calibration instead
                    %disp('Sampling frequency changed ...');
                    sampfreq(m_channel)=fread(fdat,m_length,'uint32');   % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, sampling frequency %d ...',m_channel,sampfreq(m_channel)));
                    if (chan_type(m_channel)==GL4KAUX_CHAN)
                        if (ntr(m_channel)>0)
                            trials(ntr(m_channel)).aux(m_channel).sampling_frequency = sampfreq(m_channel);
                            trials(ntr(m_channel)).aux(m_channel).time_calibration = 1e6/sampfreq(m_channel);
                        end;
                    else
                        if (ntr(m_channel)>0)
                            trials(ntr(m_channel)).channels(m_channel).sampling_frequency=sampfreq(m_channel);
                            trials(ntr(m_channel)).channels(m_channel).time_calibration=1e6/sampfreq(m_channel);
                        end;
                    end;
                case MSG_REWARD
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %disp('Reward Event record ...');   % APM
                    disp(sprintf('External Event Input OFF at %d',m_data));
                    if (ntr(m_channel)>0)
                        nevt1(m_channel)=nevt1(m_channel)+1;
                        %trials(ntr(m_channel)).channels(m_channel).events1(nevt1(m_channel))=m_data; % For APM-02 store event timestamp in the spike channel structure
                        trials(ntr(m_channel)).aux(m_channel).events1(nevt1(m_channel))=m_data; % For Guideline, store event timestamp in the auxiliary channel structure
                    end;
                case MSG_LFP_TIME_CALIBRATION
                    %disp('LFP Time Calibration record ...');
                    lfptcal(m_channel)=fread(fdat,m_length,'float32');
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).lfp_time_calibration=lfptcal(m_channel);
                    end;
                case MSG_DIGINPUT_CHANGED
                    m_data=fread(fdat,m_length,'uint32');
                    disp(sprintf(' Digital event on channel %d, timestamp %d, bits %s, bits changed %s',m_channel,m_data(1),dec2bin(m_data(2),16),dec2bin(m_data(3),16)));
                case MSG_STROBED_DATA
                    m_data=fread(fdat,m_length,'uint32');
                    disp(sprintf(' Strobed digital data on channel %d, timestamp %d, binary data %s',m_channel,m_data(1),dec2bin(m_data(2),16)));
                case MSG_SPIKE_LEN
                    spike_len(m_channel)=fread(fdat,m_length,'uint32');     % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, spike length %d ...',m_channel,spike_len(m_channel)));
                case MSG_AVG_SAMPLES
                    m_data=fread(fdat,m_length,'uint32');
                    disp(sprintf('Channel %d, chart recorder sampling frequency is %.2f ...',m_channel,1e6/(tcal(m_channel)*m_data(1))));
                case MSG_DRIVE_DEPTH
                    m_data=fread(fdat,m_length,'uint32');
                    disp(sprintf(' Drive %d depth %d ...',m_channel,m_data(1)));
                case MSG_LINE_NOISE_FILTER
                    LNF(m_channel)=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, line noise filter %d ...',m_channel,LNF(m_channel)));
%                     if (ntr(m_channel)>0)
%                         trials(ntr(m_channel)).channels(m_channel).LNF=LNF(m_channel);
%                     end;
                case MSG_HIGH_PASS_FREQ
                    high_pass_freq(m_channel)=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, high-pass filter frequency %d ...',m_channel,high_pass_freq(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).high_pass_freq=high_pass_freq(m_channel);
                    end;
                case MSG_LOW_PASS_FREQ
                    low_pass_freq(m_channel)=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    disp(sprintf('Channel %d, low-pass filter frequency %d ...',m_channel,low_pass_freq(m_channel)));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).low_pass_freq=low_pass_freq(m_channel);
                    end;
                case MSG_CHANNEL_NAME
                    m_data=fread(fdat,m_length*4,'char');    % Read message data, a string in this case
                    channel_name=char(m_data.');
                    disp(sprintf(' Channel name : %s ...',channel_name));
                case MSG_TIMESTAMP
                    %disp(sprintf(' Channel %d, %d timestamps ...',m_channel,m_length/2));
                    m_data=double(fread(fdat,m_length,'uint32'));    % Read message data, a 32-bit integer
                    if (ntr(m_channel)>0)
                        for i=1:(m_length/2)
                            %                        disp(sprintf('  %d : %x',i,m_data(2*(i-1)+1)));
                            nts(m_channel)=nts(m_channel)+1;
                            ts=(m_data(2*(i-1)+2));
                            if (ts<last_ts(m_channel))
                                disp(sprintf('Time reversal in timestamps, last value %d, new value %d !',last_ts(m_channel),ts));
                                % Assume an SOT code is missing
                                nts(m_channel)=0;
                                nwf(m_channel)=0;
                                last_ts(m_channel)=0;
                                nevt(m_channel)=0;      % Reset the event count
                                nevt1(m_channel)=0;      % Reset the event count
                                if (split)
                                    ntr(m_channel)=ntr(m_channel)+1;
                                end;
                                % Set current trial parameters that are available so far ...
                                if (~isempty(current_filename))
                                    trials(ntr(m_channel)).filename=current_filename;   % The current filename message is sent BEFORE the start-of-trial message
                                end;
                                trials(ntr(m_channel)).channels(m_channel).trig_level=trig_level(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).trig_time=trig_time(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).win_time=win_time(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).win_low=win_low(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).win_high=win_high(m_channel);
                                %trials(ntr(m_channel)).channels(m_channel).sampling_frequency=sampfreq(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).time_calibration=tcal(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).lfp_time_calibration=lfptcal(m_channel);
                                trials(ntr(m_channel)).channels(m_channel).voltage_calibration=vcal(m_channel);
                                for i=1:nunits
                                    trials(ntr(m_channel)).channels(m_channel).template(i).unit=i;
                                    try
                                        trials(ntr(m_channel)).channels(m_channel).template(i).waveform=templates(m_channel).unit(i).waveform;
                                    catch
                                        trials(ntr(m_channel)).channels(m_channel).template(i).waveform=[];
                                    end;
                                end;
                            end;
                            trials(ntr(m_channel)).channels(m_channel).timestamps(nts(m_channel))=ts;
                            last_ts(m_channel)=ts;
                        end;
                    end;
                case {MSG_NACC_WAVEFORM,MSG_ACC_WAVEFORM,MSG_NOTRIAL_NACC_WAVEFORM,MSG_NOTRIAL_ACC_WAVEFORM}
                    % if unit=0 but spike is accepted, that means we are
                    % using the window discriminator instead of template
                    % matching
                    if ((~m_unit) & ((m_code==MSG_ACC_WAVEFORM) | (m_code==MSG_NOTRIAL_ACC_WAVEFORM)))
                        m_unit=1;
                    end;
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    if (ntr(m_channel)>0)
                        nwf(m_channel)=nwf(m_channel)+1;
                    %                     disp(sprintf('%d: %s %d',nwf(m_channel), current_filename, m_length));
                        trials(ntr(m_channel)).channels(m_channel).spikes(nwf(m_channel)).unit=m_unit;
                        trials(ntr(m_channel)).channels(m_channel).spikes(nwf(m_channel)).timestamp=m_data(1);
                        trials(ntr(m_channel)).channels(m_channel).spikes(nwf(m_channel)).waveform=m_data(2:m_length);
                    end;
                    %disp(sprintf('Channel %d, waveform at %d, %d samples ...',m_channel,m_data(1),m_length-1));
                case MSG_CONTINUOUS
                    %disp(sprintf(' Continuous data packet, %d samples ...',m_length-1));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %disp(sprintf('Channel %d, continuous data packet, length %d, timestamp %d',m_channel,m_length,m_data(1)));
                    if ((ntr(m_channel)>0) && (m_data(1)>last_cts(m_channel)) && (last_cts(m_channel)>0))
                        disp(sprintf(' %d samples dropped.',m_data(1)-last_cts(m_channel)));
                        % Pad the data vector with zeros, in the same number as the missing samples
                        trials(ntr(m_channel)).channels(m_channel).continuous_samples = trials(ntr(m_channel)).channels(m_channel).continuous_samples + m_data(1)-last_cts(m_channel);
                        pause(2);
                    end;
                    last_cts(m_channel)=m_data(1)+m_length-1;   % This is in fact pointing to the expected timestamp of the next packet
                    if (ntr(m_channel)>0)
                        cs=trials(ntr(m_channel)).channels(m_channel).continuous_samples;
                        if (cs==0)
                            trials(ntr(m_channel)).channels(m_channel).start_continuous = m_data(1);    % Timestamp (in ADC samples) where recording continuous data has started
                        end;
                        if (substract_noise && ( trials(ntr(m_channel)).channels(m_channel).continuous_packets>min_cont_pkts))
                            %size(double(m_data))
                            %size(trials(ntr(m_channel)).channels(m_channel).noise_replica(1:(m_length-1)))
                            %size(trials(ntr(m_channel)).channels(m_channel).continuous)
                            trials(ntr(m_channel)).channels(m_channel).continuous(cs+(1:(m_length-1))) = double(m_data(2:end))- trials(ntr(m_channel)).channels(m_channel).noise_replica(1:(m_length-1));
                        else
                            trials(ntr(m_channel)).channels(m_channel).continuous(cs+(1:(m_length-1))) = double(m_data(2:end));
                        end;
                        trials(ntr(m_channel)).channels(m_channel).continuous_samples = trials(ntr(m_channel)).channels(m_channel).continuous_samples + m_length-1;

%                         try
%                             if (substract_noise && ( trials(ntr(m_channel)).channels(m_channel).continuous_packets>min_cont_pkts))
%                                 trials(ntr(m_channel)).channels(m_channel).continuous = cat(1,trials(ntr(m_channel)).channels(m_channel).continuous,double(m_data(2:end))-noise_replica(1:(m_length-1),m_channel));
%                             else
%                                 trials(ntr(m_channel)).channels(m_channel).continuous = cat(1,trials(ntr(m_channel)).channels(m_channel).continuous,m_data(2:end));
%                             end;
%                             % Can't assign the start_continuous only in the catch block, since the continuous field is present (though it may be empty) first time a continuous data packet is received on any channel
%                             if (isempty(trials(ntr(m_channel)).channels(m_channel).start_continuous))
%                                 trials(ntr(m_channel)).channels(m_channel).start_continuous = m_data(1);    % Timestamp (in ADC samples) where recording continuous data has started
%                             end;
%                         catch
%                             trials(ntr(m_channel)).channels(m_channel).continuous = m_data(2:m_length);
%                             trials(ntr(m_channel)).channels(m_channel).start_continuous = m_data(1);    % Timestamp (in ADC samples) where recording continuous data has started
%                         end;
                    end;
                case MSG_TEMPLATE
                    disp(sprintf(' Template waveform, unit %d, %d samples ...',m_unit,m_length-1));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    templates(m_channel).unit(m_unit).waveform=m_data(2:m_length);
                    %                     disp(sprintf('%d: %s %d',nwf(m_channel), current_filename, m_length));
                    if (ntr(m_channel)>0)
                        trials(ntr(m_channel)).channels(m_channel).template(m_unit).waveform=m_data(2:m_length);
                    end;
                case MSG_LFP
                    %disp(sprintf(' LFP data packet, %d samples ...',m_length-1));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %disp(sprintf(' LFP data packet, %d samples at %d samples ...',m_length-1,m_data(1)));
                    if (ntr(m_channel)>0)
                        try
                            trials(ntr(m_channel)).channels(m_channel).LFP = cat(1,trials(ntr(m_channel)).channels(m_channel).LFP,m_data(2:end));
                        catch
                            trials(ntr(m_channel)).channels(m_channel).LFP = m_data(2:m_length);
                        end;
                    end;
                case MSG_SLOW_DATA   % Data used for chart recorder display mode in guideline
                    %disp(sprintf(' Slow data packet, %d samples ...',m_length-1));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %disp(sprintf('Channel %d, slow data packet at %d, %d samples ...',m_channel,m_data(1),m_length-1));
                    if (ntr(m_channel)>0)
                        chan_type(m_channel)=GL4KSP_CHAN;
                        trials(ntr(m_channel)).channels(m_channel).type=chan_type(m_channel);
                        chart_data=m_data(2:end);
                        pkt_length = m_length-1;
                        try
                            trials(ntr(m_channel)).channels(m_channel).chart.min = cat(1,trials(ntr(m_channel)).channels(m_channel).chart.min,chart_data(1:2:(pkt_length-1)));
                            trials(ntr(m_channel)).channels(m_channel).chart.max = cat(1,trials(ntr(m_channel)).channels(m_channel).chart.max,chart_data(2:2:pkt_length));
                            trials(ntr(m_channel)).channels(m_channel).chart_time_calibration = trials(ntr(m_channel)).channels(m_channel).time_calibration * (2*(m_data(1)-last_chts(m_channel))/pkt_length); % The spike time calibration times the number of averaged samples
                        catch
                            trials(ntr(m_channel)).channels(m_channel).start_chart = m_data(1);    % The timestamp (in ADC samples) where chart data recording has started
                            trials(ntr(m_channel)).channels(m_channel).chart.min = chart_data(1:2:(pkt_length-1));
                            trials(ntr(m_channel)).channels(m_channel).chart.max = chart_data(2:2:pkt_length);
                        end;
                    end;
                    last_chts(m_channel)=m_data(1);
                case MSG_AUX_DATA
                    %disp(sprintf('Auxiliary I/O data on channel %d',m_channel));
                    m_data=fread(fdat,m_length,'int32');    % Read message data, a 32-bit integer
                    %m_data
                    %dec2hex(m_data)
                    %disp(sprintf('Aux data on channel %d, length %d, active channels 0x%s, timestamp %d ',m_channel,m_length,dec2hex(m_unit),m_data(1)));
                    if (ntr(m_channel)>0)
                        chan_type(m_channel)=GL4KAUX_CHAN;
                        aux_inputs = m_unit;
                        nauxin = 0;
                        while (aux_inputs>0)
                            if (bitand(aux_inputs,1)==1)
                                nauxin=nauxin+1;
                            end;
                            aux_inputs = bitshift(aux_inputs,-1);
                        end;
                        aux_data=m_data(2:end);
                        pkt_length = m_length-1;
                        for i=1:nauxin
                            try
                                trials(ntr(m_channel)).aux(m_channel).input(i).continuous = cat(1,trials(ntr(m_channel)).aux(m_channel).input(i).continuous,aux_data(i:nauxin:(pkt_length-nauxin+i)));
                                %length(trials(ntr(m_channel)).aux(m_channel).input(i).continuous)
                            catch
                                trials(ntr(m_channel)).aux(m_channel).input(i).start_continuous = m_data(1);    % The timestamp (in ADC samples) where continuous data recording has started
                                trials(ntr(m_channel)).aux(m_channel).input(i).continuous = aux_data(i:nauxin:(pkt_length-nauxin+i));
                            end;
                        end;
                    end;
                case MSG_TCPIP_USER1
                    m_data=fread(fdat,m_length*4,'char');    % Read message data, a string in this case
                    current_filename=char(m_data.');
                    disp(sprintf(' Trial Filename : %s ...',current_filename));
                otherwise
                    %length(m_data)
                    %dec2hex(m_data) % Previous data
                    disp(sprintf('Unknown message code %d on channel %d, length %d at file position %d (%d DWORDS).',m_code,m_channel,m_length,ftell(fdat),ftell(fdat)/4));
                    m_data=fread(fdat,m_length,'uint32');    % Read message data, a 32-bit integer
                    %dec2hex(m_code)
                    %[dec2hex(m_unit) dec2hex(m_channel)]
                    %dec2hex(m_length)
                    %dec2hex(m_data(1:1000))
                    %pause;
                    
            end;
        end;
    end;
    fclose(fdat);
    % The self-interference (or from channels with correlated activity) noise has been substracted earlier
    % Now substract noise that occurs cross-channels, for instance from channels having different sampling rate
    if (substract_noise)
        for i=1:nchan
            if (chan_type(i)==GL4KSP_CHAN)
                for j=1:ntr(i)
                    if (length(trials(j).channels)>=i)
                        %disp(sprintf('Channel %d, trial %d of %d',i,j,ntr(i)));
                        if (trials(j).channels(i).continuous_samples>0)
                            if (isempty(fixed_noise_freq))
                                pkt_intvls = [];
                            else
                                pkt_intvls = round(sampfreq(i)./fixed_noise_freq);
                            end;
                            for k=1:nchan
                                if (chan_type(k)==GL4KSP_CHAN)
                                    if (length(trials(j).channels)>=k)
                                        % Check for data packet intervals OTHER than the one used for the current channel
                                        if (trials(j).channels(k).continuous_samples>0)
                                            pkt_intvl = round(pkt_len(k)*sampfreq(i)/sampfreq(k));
                                            if (pkt_intvl~=pkt_len(i))
                                                pkt_intvls = [pkt_intvls pkt_intvl];
                                                disp(sprintf('Channel %d, trial %d, artifact every %d samples from channel %d',i,j,pkt_intvl,k));
                                            end;
                                        end;
                                    end;
                                end;
                                if (chan_type(k)==GL4KAUX_CHAN)
                                    if (isfield(trials(j),'aux') && (length(trials(j).aux)>=k))
                                        % Check for noise introduced by the aux channels as well
                                        m=1;    % Aux input 1 - all inputs share the same sampling frequency
                                        naux_inputs=2; % Assume two inputs
                                        if (length(trials(j).aux(k).input(m).continuous)>0)
                                            pkt_intvl = round(pkt_len(k)*sampfreq(i)/sampfreq(k))*naux_inputs;
                                            if (pkt_intvl~=pkt_len(i))
                                                pkt_intvls = [pkt_intvls pkt_intvl];
                                                disp(sprintf('Channel %d, trial %d, artifact every %d samples from aux channel %d',i,j,pkt_intvl,k));
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                            if (~isempty(pkt_intvls))
                                pkt_intvls(find(pkt_intvls < 2)) = []; % Replica must be at least 2 samples long
                                pkt_intvls=sort(unique(pkt_intvls));
                                %disp(sprintf('  Channel %d, trial %d, artifacts every:',i,j));
                                %disp(sprintf('\t%d',pkt_intvls));
                                for n=1:length(pkt_intvls)
                                    num_packets=floor(trials(j).channels(i).continuous_samples/pkt_intvls(n));
                                    if (num_packets>min_cont_pkts)
                                        % Truncate the wave to a multiple of the packet interval
                                        wav=trials(j).channels(i).continuous(1:(num_packets*pkt_intvls(n)));
                                        %nrep=zeros(pkt_intvls,1);
                                        nwav=reshape(wav,[pkt_intvls(n) num_packets]);
                                        nwav=mean(nwav,2);  % This is the noise replica
                                        %figure(77);clf;plot(nwav);pause;
                                        ix=1:pkt_intvls(n);
                                        for i1=0:(num_packets-1)
                                           trials(j).channels(i).continuous(ix)=trials(j).channels(i).continuous(ix)-nwav;
                                           ix=ix+pkt_intvls(n);
                                        end;
                                        % Check if there's a last (truncated) packet
                                        ix=(num_packets*pkt_intvls(n)+1):length(trials(j).channels(i).continuous);
                                        if (~isempty(ix))
                                           trials(j).channels(i).continuous(ix)=trials(j).channels(i).continuous(ix)-nwav(1:length(ix));
                                        end;
                                    end;
                                end;
                            end;
                        end;
%                         if (isfield(trials(j).channels(i),'continuous_packets'))
%                             if (trials(j).channels(i).continuous_packets>min_cont_pkts)
%                                 trials(j).channels(i).noise_replica = trials(j).channels(i).noise_replica/trials(j).channels(i).continuous_packets;
%                             end;
%                         end;
                    end;
                end;
            end;
        end;
    end;
    
else
    disp(sprintf('Could not open %s ...',filename));
end;

apmdata=trials;