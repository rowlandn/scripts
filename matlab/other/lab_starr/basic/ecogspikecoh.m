
% This program calculates and plots the coherence between a spike train
% (timestamps of spikes from single unit recording) and 5 simultaneously
% recorded bipolar ecog recordings. optional:It also plots the power spectral
% density of each ecog channel and of the spike train using the "spike
% shuffling method".  It is designed for data input generated from one type
% of movement only (all rest or all active, not alternating rest/active)
% Inputs: 
%       spk.t = vector of times of spikes (in seconds, which must be imported
%           from a nexfile that was created by plexon sorting of single
%           unit data(note that maybe plexon can save spike times in a mat
%           file or somehting very easy to open in matlab)
%       ecog_lfp_raw_data = matrix of raw ecog data, with first six columns
%           representing ecog recordings, generated using apmconv7_coh

% remontage the ecog channels to bipolar recordings

%load a nexfile containing spiketimes - note rob's program runspikeoscil
%has a lengthy section that does this, or check plexon what file formats it
%can save things as

%% Define Variables
Fs=1000;           
NFFT = 512;
window = 512; %epoch data will be divided into chunks of size equal to window for coherence calculation
noverlap = 256;
global N_SHUF;
N_SHUF = 100;
%% Import ecog data
% import filename and pathname of mat file containing ecog data created by
% apmconv7_coh
[fn pn] = uigetfile('*.mat','Select .mat containing _ecog data'); 
cd(pn);
load([pn fn]);
% remove '_ecog_lfp.mat' ending from filename
fn = strrep(fn,'_ecog.mat','');

%% Remontage Ecog Data
%such that contact 1 is referencing contact 2, 2 references 3, etc

% Recog = zeros(size(ecog_lfp_raw_data(:,end-2);
ecog1v2 = ecog.contact_pair(1,1).raw_ecog_signal - ecog.contact_pair(1,2).raw_ecog_signal;
ecog2v3 = ecog.contact_pair(1,2).raw_ecog_signal - ecog.contact_pair(1,3).raw_ecog_signal;
ecog3v4 = ecog.contact_pair(1,3).raw_ecog_signal - ecog.contact_pair(1,4).raw_ecog_signal;
ecog4v5 = ecog.contact_pair(1,4).raw_ecog_signal - ecog.contact_pair(1,5).raw_ecog_signal;
ecog5v6 = ecog.contact_pair(1,5).raw_ecog_signal;

recog = [ecog1v2; ecog2v3; ecog3v4; ecog4v5; ecog5v6];
[num_row num_col]=size(recog);
num_contact_pair = num_row;
%% Import Spike data
%obtains number of spikes and spike timestamps from .nex file, taking
%advantage of code found in nex2isi.m

global RTPATH	% starting directory containing NEX file

[nexname RTPATH] = uigetfile('*.nex','Select spike time file (NEX)');

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

fname = strrep(nexname,'.nex','');
isi = 1000 .* diff(spk.t);
bad = find(isi < 1);
if length(bad)
    warning(['Found ' int2str(length(bad)) ' ISI''< 1 msec!!  Correcting...']);
    spk.t( bad+1 ) = [];
    spk.n = length(spk.t);
end

%% Convert spike times to delta function as per Halliday 1996 article
spk.t = round( 1000.*spk.t );
len_delt = spk.t(end);
spk.delt = zeros(1,len_delt);
spk.delt(spk.t) = 1; %array indexing to assign spk.delt =1 whenever there is a spike
% %% Downsampling the single unit
% %MER done at 20K Hz, whereas ecog at 1K Hz
% 
% spk.down = downsample(spk.delt,20);
% len_down = length(spk.down);
%%
% MAKE SURE WHEN YOU USE MSCOHERE OR PSD OR OTHER FREQUENCY DOMAIN
% TRANSFORM ON THE SPK.DELT FILE, YOU HAVE TO SUBTRACT THE MEAN VALUE OF
% THE FILE TO AVOID LOW FREQUENCY ARTIFACTS IN THE POWER SPECTRUM (OR
% COHERENCE SPECTRUM)

% Calculate the power on the spike train (this uses Rob's old method with
% psd function and spike shuffling
hann = hanning(NFFT); %we will use pwelch for psd which uses hanning windows automatically
[pow, freq] = psd( spk.delt-mean(spk.delt), NFFT, Fs, hann); %this is "linear detrending"
% [pow, freq] = psd( spk.down-mean(spk.down), NFFT, Fs, hann); %this is "linear detrending"
spk.pow = pow';
spk.freq = freq';
len_spect = length(freq);

% Calculate psd's from globaly shuffled spike trains with the same ISI
fprintf('Shuffling');
clk_tk = round(N_SHUF/10);
isi = diff(spk.t);
isi = isi(isi>0); %make sure isi's are good
pow_shuf=zeros( N_SHUF, len_spect);
for i=1:N_SHUF
    r=randperm(length(isi));
    rand_isi(r)=isi;
    y=cumsum(rand_isi);
    rand_delt = zeros(1,len_delt);
    rand_delt([1 y+1])=1; 
    [pow_rand,freq_rand] = psd( rand_delt-mean(rand_delt), NFFT, Fs, hann);
    pow_shuf(i,:)=pow_rand';
	if ~mod(i,clk_tk);	fprintf('.');	end
end
% for i=1:N_SHUF
%     r=randperm(length(isi));
%     rand_isi(r)=isi;
%     y=cumsum(rand_isi);
%     rand_down = zeros(1,len_down);
%     rand_down([1 y+1])=1; 
%     [pow_rand,freq_rand] = psd( rand_down-mean(rand_down), NFFT, Fs, hann);
%     pow_shuf(i,:)=pow_rand';
% 	if ~mod(i,clk_tk);	fprintf('.');	end
% end
spk.pow_rand = mean(pow_shuf);	% Mean pow of randomized ISI's
fprintf('\n');
%%
% Now compute ISI-compensated normalized spectrum
spk.pow_comp = spk.pow ./ spk.pow_rand;

%% Calculate PSD of spike data
spk.delt_minusmean = spk.delt-mean(spk.delt); %removing low frequency artifact
[PSDspk,F] = pwelch(spk.delt_minusmean,512,256,512,1000);
hf = figure('Name','PSD spike data');
    subplot(2,1,1)
    plot(F,PSDspk)
    xlim([0 150]);
        if i==1;
            str=[fn sprintf('\n') 'PSD spike data']; % allows title to have file name
        else
            str=['PSD spike data'];
        end
        title(str);
     subplot (2,1,2)
     plot(F,PSDspk)
     xlim([0 50]);
        if i==1;
            str=[fn sprintf('\n') 'PSD spike data']; % allows title to have file name
        else
            str=['PSD spike data'];
        end
        title(str);
%% Calculate PSD of ecog data
%Should be identical to Sho's code using spectrogram function, no? Why use
%spectrogram here and psd for spike data? Are they equivalent?

[PSD1,F] = pwelch(recog(1,:),512,256,512,1000);
[PSD2,F] = pwelch(recog(2,:),512,256,512,1000);
[PSD3,F] = pwelch(recog(3,:),512,256,512,1000);
[PSD4,F] = pwelch(recog(4,:),512,256,512,1000);
[PSD5,F] = pwelch(recog(5,:),512,256,512,1000);

PSDecog = [PSD1 PSD2 PSD3 PSD4 PSD5];

hf = figure('Name','PSD ecog');
for i=1:5
        subplot(2,5,i)
        plot(F,PSDecog(:,i))
        xlim([0 150]);
        if i==1;
            str=[fn sprintf('\n') 'PSDecog' num2str(i)]; % allows title to have file name
        else
            str=['PSDecog' num2str(i)];
        end
        title(str);
end
hold on
for i=1:5
        subplot(2,5,i+5)
        plot(F,PSDecog(:,i))
        xlim([0 50]);
        if i==1;
            str=[fn sprintf('\n') 'PSDecog' num2str(i)]; % allows title to have file name
        else
            str=['PSDecog' num2str(i) ];
        end
        title(str);
end
%% Calculate Coherence between ecog and spike train. 
% The two signals have different sampling rates. Need to down-sample spike data?
Coh = zeros(NFFT/2+1, num_contact_pair);

spk.delt_minusmean = spk.delt-mean(spk.delt); %removing low frequency artifact
% spk.down_minusmean = spk.down-mean(spk.down); %removing low frequency artifact
len_ecog = length(recog);

%MSCOHERE requires two signals of equal length
if len_ecog>len_delt

    for i = 1:num_contact_pair;
%         
        inp1 = spk.delt_minusmean;
        inp2 = recog(i,1:len_delt); 

        [Cxy,F] = mscohere(inp1,inp2,window,noverlap,NFFT,Fs);
        Coh(:,i) = Cxy;
    end

elseif len_delt>len_ecog

    for i = 1:num_contact_pair;
        inp1 = spk.delt_minusmean(1:len_ecog);
        inp2 = recog(i,:); 

        [Cxy,F] = mscohere(inp1,inp2,window,noverlap,NFFT,Fs);
        Coh(:,i) = Cxy;
    end
    
else

    for i = 1:num_contact_pair;
        inp1 = spk.delt_minusmean;
        inp2 = recog(i,:); 

        [Cxy,F] = mscohere(inp1,inp2,window,noverlap,NFFT,Fs);
        Coh(:,i) = Cxy;
    end
end
%% Calculate inverse hyperbolic tangent for transformed coherence 
% Y = atanh(X) returns the inverse hyperbolic tangent for each element of X.

trans_coh = atanh(sqrt(Coh));

%% Plotting Coherence
% hf = figure('Name',(fn));
% 
% for i = 1:num_contact_pair
%     h = subplot(5,1,i);
%     plot(F, trans_coh(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherencet','color','b','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 30 0 0.5]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
% 
%    if i==1
%        ylabel('transformed Coherence')
%        xlabel('F')
%    end
% 
% end
%% Plotting coherence
% hf = figure('Name',(fn));
% 
% for i = 1:num_contact_pair
%     h = subplot(2,5,i);
%     plot(F, Coh(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherencet','color','b','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 30 0 0.5]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
%    
%     if i==1
%        ylabel('Coherence')
%        xlabel('F')
%     end
% 
%     h = subplot(2,5,i+5);
%     plot(F, trans_coh(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherencet','color','b','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 30 0 0.5]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
%    
%    if i==1
%        ylabel('transformed Coherence')
%        xlabel('F')
%    end
% end

% createfigure(F,Coh(:,1),trans_coh(:,1),Coh(:,2),trans_coh(:,2),Coh(:,3),trans_coh(:,3),Coh(:,4),trans_coh(:,4),Coh(:,5),trans_coh(:,5))
% %CREATEFIGURE(F1,MEAN COHERENCET1,MEAN COHERENCET2,MEAN COHERENCET3,MEAN COHERENCET4,MEAN COHERENCET5,MEAN COHERENCET6,MEAN COHERENCET7,MEAN COHERENCET8,MEAN COHERENCET9,MEAN COHERENCET10)
% %  F1:  vector of x data
%  MEANCOHERENCET1=Coh(:,1)
%  MEANCOHERENCET2=trans_coh(:,1)
%  MEANCOHERENCET3=Coh(:,2)
%  MEANCOHERENCET4=trans_coh(:,2)
%  MEANCOHERENCET5=Coh(:,3)
%  MEANCOHERENCET6=trans_coh(:,3)
%  MEANCOHERENCET7=Coh(:,4)
%  MEANCOHERENCET8=trans_coh(:,4)
%  MEANCOHERENCET9=Coh(:,5)
%  MEANCOHERENCET10=trans_coh(:,5)
% 
% %  Auto-generated by MATLAB on 24-Nov-2008 11:03:31
% 
% Create figure
figure1 = figure('Name',(fn));

% Create axes
axes1 = axes('Parent',figure1,'Position',[0.06816 0.7727 0.3963 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])

box('on');
hold('all');

% Create plot
plot(F,Coh(:,1),'DisplayName','Rest','Parent',axes1,...
    'LineWidth',1.5);

% Create ylabel
ylabel('Coherence');

% Create xlabel
xlabel('F');

% Create axes
axes2 = axes('Parent',figure1,'Position',[0.5341 0.7727 0.3969 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,trans_coh(:,1),'DisplayName','Rest','Parent',axes2,...
    'LineWidth',1.5);

% Create ylabel
ylabel('transformed Coherence');

% Create xlabel
xlabel('F');

% Create axes
axes3 = axes('Parent',figure1,'Position',[0.06816 0.6034 0.4 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,Coh(:,2),'DisplayName','Rest','Parent',axes3,...
    'LineWidth',1.5);

% Create axes
axes4 = axes('Parent',figure1,'Position',[0.5341 0.6034 0.3969 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,trans_coh(:,2),'DisplayName','Rest','Parent',axes4,...
    'LineWidth',1.5);

% Create axes
axes5 = axes('Parent',figure1,'Position',[0.06816 0.4341 0.3969 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,Coh(:,3),'DisplayName','Rest','Parent',axes5,...
    'LineWidth',1.5);

% Create axes
axes6 = axes('Parent',figure1,'Position',[0.5341 0.4341 0.3969 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,trans_coh(:,3),'DisplayName','Rest','Parent',axes6,...
    'LineWidth',1.5);

% Create axes
axes7 = axes('Parent',figure1,'Position',[0.06816 0.2648 0.3962 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,Coh(:,4),'DisplayName','Rest','Parent',axes7,...
    'LineWidth',1.5);

% Create axes
axes8 = axes('Parent',figure1,'Position',[0.5341 0.2648 0.3969 0.1209]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,trans_coh(:,4),'DisplayName','Rest','Parent',axes8,...
    'LineWidth',1.5);

% Create axes
axes9 = axes('Parent',figure1,'Position',[0.06944 0.09553 0.3973 0.1257]);
axis([0 150 0 0.25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,Coh(:,5),'DisplayName','Rest','Parent',axes9,...
    'LineWidth',1.5);

% Create axes
axes10 = axes('Parent',figure1,'Position',[0.5349 0.09553 0.3978 0.1247]);
axis([0 150 0 .25])
set (gca,'xtick',[0:20:150])
box('on');
hold('all');

% Create plot
plot(F,trans_coh(:,5),'DisplayName','Rest','Parent',axes10,...
    'LineWidth',1.5);