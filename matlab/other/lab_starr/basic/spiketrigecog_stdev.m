function spiketrigecog_stdev()
% spiketrigecog()
% This m-file takes single-unit spike timestamps stored in a nex, then
% runs spike-triggered signal averaging and on ecog recordings.  The
% program helps detect proper standard deviation level

%
%
% Modified version of spiketrigecog.m, SAS 6/29/2010

%% Initialize varaibles
Fs = 1000;      % GS4k: sampling rate for ecog
% Fs = 751;      % AlphaOmega: sampling rate for ecog
PRE_SPK = 0.5; % pre-spike trigger period (sec)
PST_SPK = 0.5;  % post-spike trigger period (sec)
STDEV = 3;     % standard deviation for significance
NSPK = 1000;     % number of spikes to use for signal averaging

%% spike timestamps
[nexfn, nexpn] = uigetfile('*.nex', 'Select spike time file (NEX)');
filename = strrep(nexfn,'.nex','');

if (nexfn ~= 0)
    cd(nexpn);
    [nvar, varname, types] = nex_info(nexfn);
    if nvar > 1      
        % allow user to select desired channel if there are multiple
        % channels
        varnamecell = cellstr(varname);
        varnameidx = menu(['There are > 1 channels in ' nexfn '.' sprintf('\n')...
        'Choose the channel name with desired spike timstamps.'], varnamecell);
        varname = varnamecell{varnameidx};
    end
    [spk.n, spk.t] = nex_ts(nexfn,varname);
else
    error(['I can''t find the NEX file:  ' nexfn ' in ' nexpn]);
end

% find ISI < 1 msec
isi = 1000 .* diff(spk.t);
bad = find(isi < 1);
if ~isempty(bad)
    warning(['Found ' int2str(length(bad)) ' ISI''< 1 msec!!  Correcting...']);
    spk.t( bad+1 ) = [];
    spk.n = length(spk.t);
end


%% ecog data

% load data
% [ecogfn, ecogpn] = uigetfile('*_ecog.mat', 'Select ecog file (*_ECOG.MAT)');
% cd(ecogpn);
ecogfn = strcat(filename,'_ecog.mat');
load(ecogfn);

% remontage (units are in microv)
ecog12= (ecog.contact_pair(1).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal);
ecog23= (ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal);
ecog34= (ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal);
ecog45= (ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal);
ecog56= ecog.contact_pair(5).raw_ecog_signal;

% pt Solis 3/23/10 recorded with ecog signal and common reference 6 flipped
% around.  correct for flipped signal.
% ecog12= (ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(1).raw_ecog_signal);
% ecog23= (ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal);
% ecog34= (ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal);
% ecog45= (ecog.contact_pair(5).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal);
% ecog56= -ecog.contact_pair(5).raw_ecog_signal;

% % AlphaOmega (SAS 2/9/10): there is an annoying DC offset in all ecog channels.
% % Eliminate offset from remontaged signals by subtracting the average value.
% ecog12=ecog12-mean(ecog12);
% ecog23=ecog23-mean(ecog23);
% ecog34=ecog34-mean(ecog34);
% ecog45=ecog45-mean(ecog45);
% ecog56=ecog56-mean(ecog56);
% 
% % AlphaOmega (SAS 2/9/10): remontaged channels have significantly lower signal amplitude
% % than ecog56, making it difficult to visualize STA waveforms on the plot.
% % Multiply remontaged channels and computationally increase their signal
% % amplitude.
% ecog12=ecog12*10;
% ecog23=ecog23*10;
% ecog34=ecog34*10;
% ecog45=ecog45*10;

ecogall = [ecog12; ecog23; ecog34; ecog45; ecog56];


%% spike-triggered average of ECog signal

spk_t = spk.t;


[ecg_n, ecg_l] = size(ecogall); % find number and length of ecog channels
% discard spk ts at beginning of recording with inadequate pre-spike period
spk_t = spk_t(spk_t > PRE_SPK); 
% discard spk ts at end of recording with inadequate post-spike period
spk_t = spk_t((spk_t+PST_SPK)*Fs < ecg_l); 
spk_n = length(spk_t);

% optional: use predetermined number of spikes
if true
    if spk_n > NSPK
        spk_t(NSPK+1:end) = []; % only include the first NSPK spikes
        spk_n = NSPK;
    else
        error([nexfn ' does not exceed the # spike cutoff of' num2str(NSPK)...
            '.  Lower minimum cutoff for this unit.']);
    end
end

% create time vector
t = 1000*(-PRE_SPK:1/Fs:PST_SPK); % multiplied by 1000 to change to msec scale

% parse ecog centered at spike ts from each contact pair
E = zeros(spk_n,length(t),ecg_n); %3D matrix for storing spike triggered ecog data
for i = 1:spk_n
    tmp1 = int32((spk_t(i)-PRE_SPK) * Fs);  % int32 used to keep index in integer format
    tmp2 = int32((spk_t(i)+PST_SPK) * Fs);
    % since int32 rounds to the next closest integer, there may be some
    % cases where the parsing indeces are off by 1.  fix them on the fly.
    d = tmp2 - tmp1;
    % lengt(t)-d must equal 1 for the parsed ecog data to fit E 3D matrix
    if length(t)-d == 0
        tmp2 = tmp2-1;
    elseif length(t)-d == 2
        tmp2 = tmp2+1;
    end
    for k = 1:ecg_n
        E(i,:,k) = ecogall(k,tmp1:tmp2);
    end
end

STA = mean(E,1);  % spike-triggered ecog 

% create randomly triggered ecog 
tend = spk_t(end)-PST_SPK; % random timestamps are bound in time between 0 sec and time of last spike minus post-spike period
rnd_t = ceil(1000*(tend*rand(spk_n,1))); % random time (msec)
rnd_t = rnd_t/1000+PRE_SPK; %random time (sec), plus pre-spike period

% parse ecog centered at spike ts from each contact pair
Er = zeros(spk_n,length(t),ecg_n); %3D matrix for storing randomly triggered ecog data
for i = 1:spk_n
    tmp1 = int32((rnd_t(i)-PRE_SPK) * Fs);  % int32 used to keep index in integer format
    tmp2 = int32((rnd_t(i)+PST_SPK) * Fs);
    % since int32 rounds to the next closest integer, there may be some
    % cases where the parsing indeces are off by 1.  fix them on the fly.
    d = tmp2 - tmp1;
    % lengt(t)-d must equal 1 for the parsed ecog data to fit E 3D matrix
    if length(t)-d == 0
        tmp2 = tmp2-1;
    elseif length(t)-d == 2
        tmp2 = tmp2+1;
    end
    for k = 1:ecg_n
        Er(i,:,k) = ecogall(k,tmp1:tmp2);
    end
end

RTA = mean(Er,1); % randomly triggered ecog average
RTA_std = std(RTA); % standard deviation of RTA

% %% plot figure
% hf1 = figure;
% hold on
% % figure 1: non-normalized, same plot
% 
% create a vector of ones used later for SD plotting
vones = ones(1,length(t));

% specify where along the time axis to place ecog pair label
t_text = -PRE_SPK*1e3*0.9;
% 
% % find min and max of STA
% vmin = min(min(min(STA)));
% vmax = max(max(max(STA)));
% 
% % calculate stacking constant used for plotting
% if (vmax-vmin) > 2*STDEV*max(RTA_std)
%     C_stk = vmax-vmin;
% else
%     %The 'else' case prevents stdev lines of neighboring ecog pairs from
%     %overlapping.
%     C_stk = 2*STDEV*max(RTA_std); 
% end
% 
% for k = 1:ecg_n
%     stk = (ecg_n - k) * C_stk;
%     z = STA(1,:,k) + stk; 
%     % create lines for STA average, upper and lower thresholds for plotting
%     avg = (STA_mean(1,1,k)*vones) + RTA(1,:,k) + stk;
%     thu = (STA_mean(1,1,k)*vones) + (STDEV * RTA_std(1,1,k)) + stk;
%     thl = (STA_mean(1,1,k)*vones) - (STDEV * RTA_std(1,1,k)) + stk;
%     plot(t,z,'LineWidth',2); % plot STA
% %     plot(t,avg,'-k'); %plot randomly triggered average
%     plot(t,thu,'-.r',t,thl,'-.r'); % plot lower and upper thresholds
%     if k==M1_ch
%         text(t_text,thu(1),...
%             ['ecog' num2str(k) num2str(k+1) ',M1'],...
%             'VerticalAlignment','Top','FontWeight','bold'); % label M1 ecog pair
%     else
%         text(t_text,thu(1),...
%             ['ecog' num2str(k) num2str(k+1)],...
%             'VerticalAlignment','Top'); % label all other ecog pairs
%     end
% end
% xlabel('Time (millisec)');
% ylabel('Non-normalized STA (microV)');
% title([filename sprintf('\n')...
%     '# spks=' num2str(spk_n) ', sig thresh=' num2str(STDEV) 'SD'],...
%     'FontSize',12);
% ylm=ylim;
% plot([0 0],[ylm(1) ylm(2)],'k--');
% hold off


%% M1 STA processing

m1sta=STA(1,:,M1_ch);
m1rta=RTA(1,:,M1_ch);

% subtract DC offset
m1sta=m1sta-mean(m1sta);
m1rta=m1rta-mean(m1rta);

% find max values for STA and 
[vmax vi]=max(m1sta);
vt=t(vi);
[vrmax vri]=max(m1rta);
vrt=t(vri);

% More plots
figure;
plot(t,m1sta),hold on
plot(t,m1rta,'k','LineWidth',2)
thu = vones*(STDEV * RTA_std(1,1,M1_ch));
thl = -vones*(STDEV * RTA_std(1,1,M1_ch));
plot(t,thu,'-.r',t,thl,'-.r'); % plot lower and upper thresholds
plot(vt,vmax,'o');
plot(vrt,vrmax,'ok');
ylm=ylim;
plot([0 0],[ylm(1) ylm(2)],'k--'); % plot vertical bar for t=0
xlabel('Time (millisec)');
ylabel('Voltage (microV)');
title([filename sprintf('\n')...
    '# spks=' num2str(spk_n) ', sig thresh=' num2str(STDEV) 'SD'],...
    'FontSize',12);
text(t_text,thu(1),...
    ['ecog' num2str(M1_ch) num2str(M1_ch+1) ', M1'],...
    'VerticalAlignment','bottom',...
    'FontWeight','bold'); % label all other ecog pairs
text(250,thl(1),...
    ['STA peak @ t=' num2str(vt) 'ms, max=' num2str(vmax) sprintf('\n')...
    'RTA peak @ t=' num2str(vrt) 'ms, max=' num2str(vrmax) sprintf('\n')],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');
legend('M1 STA','M1 RTA');
legend('boxoff');

hold off

% %% save variables and figures
% save([filename '_sta'],'t','STA','STA_mean',...
%     'm1sta','m1stafilt','xc','env','vmax','vt','vrmax','vrt','vt0','phase');
% saveas(hf1,[filename '_sta'],'fig');
% saveas(hf2,[filename '_m1sta'],'fig');



