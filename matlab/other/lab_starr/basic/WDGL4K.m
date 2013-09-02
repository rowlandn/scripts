function tsorted = WDGL4K( t, chan, discriminator_mode, trig_level, pre_trigger, win_low, win_high, win_delay, waveform_len )
% WD - Offline window discriminator
% The discrimination parameters are extracted from the trial structure, but if
% specified as arguments, override the trial values
%
% discriminator_mode values:
%       0 - window discriminator
%       1 - level discriminator (default)
%
% AB, 2008

try
    discriminator_mode;
catch
    discriminator_mode = 1;
end;

% Window discrimination parameters
try
    trig_level;
catch
    %trig_level=-0.32;    % wavread returns normalized audio data, in the interval -1 to +1
    trig_level=t.channels(chan).trig_level;
end
try
    pre_trigger;     % in audio samples
catch
    pre_trigger=t.channels(chan).trig_time;     % in audio samples
end;
try
    trig_slope;
catch
    trig_slope=t.channels(chan).trig_slope;
end
try
    win_low;
catch
    win_low = t.channels(chan).win_low;
end;
try
    win_high;
catch
    win_high = t.channels(chan).win_high;
end;
try
    win_delay;
catch
    win_delay = t.channels(chan).win_time;     % in ADC samples
end;
try
    waveform_len;
catch
    waveform_len = 64;  % in ADC samples
end;

try
    t;
    disp('Data is already loaded into memory. Type ''clear all'' to force reloading data.');
catch
    disp('Loading sample data ...');
    t=APMReadData('GL4KData.apm');
end;

y = t.channels(chan).continuous;
sf = t.channels(chan).sampling_frequency;
minspike=1.1*min(y);
maxspike=1.1*max(y);


% Now look for spikes
ts=[];      % Initialize array of spike timestamps, for raster plot

% Create a new figure where the waveforms and window discriminator
% parameters
% figure(11);
% clf;        % Clear figure
% colordef(11,'black');    % Set white-on-black color scheme
% box on;
% hold on;
% xlabel('ADC sample');
% ylabel('Normalized amplitude');
% axis([1 waveform_len minspike maxspike]);
% % Plot the discriminator parameters (trig, window)
% plot([1 waveform_len],[trig_level trig_level],'y-','LineWidth',2);
% plot([pre_trigger+win_delay+1 pre_trigger+win_delay+1],[win_low win_high],'w-','LineWidth',4);
t.channels(chan).timestamps=[];
t.channels(chan).spikes=[];
nspikes=0;
for i=(1+pre_trigger):(length(y) + pre_trigger - waveform_len)
    % Check trigger level crossing
    if (((trig_slope>0) && (y(i)<trig_level) && (y(i+1)>=trig_level)) || ((trig_slope<=0) && (y(i)>trig_level) && (y(i+1)<=trig_level)))
        acc = 0;
        if (discriminator_mode==0)
            % Check window crossing
            if ((y(i+win_delay)>win_low) && (y(i+win_delay)<win_high))
                acc=1;
            end;
        else
            % Check level crossing
            maxsp=max(y((i+1):(i+waveform_len-pre_trigger)));
            minsp=min(y((i+1):(i+waveform_len-pre_trigger)));
            %[max(y) maxsp win_low win_high]
            if (win_high>0)
                if ((maxsp>win_low) && (maxsp<win_high))
                    acc=1;
                end;
            else
                if ((minsp<win_low) && (minsp>win_high))
                    acc=1;
                end;
            end;
        end;
        if (acc)
            %plot(y((i-pre_trigger):(i-pre_trigger+waveform_len)),'g-');
            disp(sprintf(' Spike detected at %d (%.3f s) ...',i,i/sf));
            ts=[ts;i/sf];
            nspikes=nspikes+1;
            t.channels(chan).timestamps(nspikes)=t.channels(chan).start_continuous+i-1;
            t.channels(chan).spikes(nspikes).timestamp=t.channels(chan).timestamps(nspikes);
            t.channels(chan).spikes(nspikes).unit=1;
            t.channels(chan).spikes(nspikes).waveform=y((i-pre_trigger):(i-pre_trigger+waveform_len));
            i=i+waveform_len-pre_trigger;   % Skip samples to the end of the spike
        else
            % Do not save non-discriminated waveform
            %plot(y((i-pre_trigger):(i-pre_trigger+waveform_len)),'r-');
        end;
    end;
    if (mod(i,48000)==0)
        disp(sprintf(' Processing sample %d (%.1f%%) ...',i,100.0*i/length(y)));
        %drawnow;
    end;
end;
%hold off;

tsorted=t;

% % Figure 2, raster plot and histogram of the spikes
% figure(2);
% clf;        % Clear figure
% colordef(2,'black');    % Set white-on-black color scheme
% subplot(2,1,1);
% box on;
% hold on;
% title('Raster Plot');
% xlabel('Time (s)');
% axis([0 length(y)/sf 0 1]);
% for i=1:length(ts)
%     plot([ts(i) ts(i)],[0 1],'g-','LineWidth',4);
% end;
% hold off;
% 
% subplot(2,1,2);
% box on;
% hold on;
% title('Histogram');
% ylabel('Number of spikes');
% xlabel('Time (s)');
% hist(ts,200);
% hold off;
% 
% 
% 
% % Average firing rate over the selected interval
% mfr=length(ts)/(length(y)/sf);
% disp(sprintf('Mean firing rate over the selected interval is %.2f sp/sec',mfr));
% 
% % Figure 3, interspike interval
% figure(3);
% clf;        % Clear figure
% colordef(3,'black');    % Set white-on-black color scheme
% box on;
% hold on;
% title('ISI Histogram');
% ylabel('Number of spikes');
% xlabel('ISI (msec)');
% ISI=diff(ts);
% hist(1000*ISI,100); % Multiply with 1000 to get result in msec
% % mean ISI (should match 1/mfr)
% mISI=mean(ISI);
% disp(sprintf('Mean ISI is: %.2f msec',1000*mISI));  % Multiply with 1000 to get result in msec
% hold off;
