%% ecogCOH
%This program is based on recog_coh_all_data program for calculating
%coherence of ecog and lfp signals. Where recog_coh_all_data used matrix
%format, ecogCOH uses structure format and is designed to be parallel to
%the ecogPSD program. 

%INPUT: Like ecogPSD it uses *_ecog.mat data, which is the output of apmconv7. 

%OUTPUT: This program outputs transformed coherence data as *_trans_coh.mat
%files for later group coherence analysis.

%ALC 6/5/09

% ALC 6/10/09: updated to process M1-S1 coherence simultaneously. This data
% is also stored in the transcoh structure for output and later group
% analysis. 

% ALC 6/19/09: updated M1_ch references to be in line with update to
% apmconv7 code. 
%% Define Constants

OFFSET = 0.5; % epoch time offset in seconds, NOTE we may change this to zero
Fs=1000;            % samples per unit time, in this case digitization at 1000 hz rate
% Fs = 1.502403855323792e3; %Sampling rate of alphaOmega system
nfft = 512; %FFT length that determines the frequencies at which the coherence is estimated. For real x and y,the length of Cxy is (nfft/2+1) if nfft is even or (nfft+1)/2 if nfft is odd. 
window = 512; %epoch data will be divided into chunks of size equal to window for coherence calculation
noverlap = 256; % how much overlap there is between windows
% nfft and window must be equal; noverlap is a fraction of those (50%,75%, etc.)
EPOCH_LEN = 2; % in seconds

%% Import/Prepare file
% import filename and pathname of .mat file containing ecog data created by
% apmconv7
[fn pn] = uigetfile('*_ecog.mat','Select .mat containing _ecog data');
cd(pn);
load([pn fn]);
% remove '_ecog.mat' ending from filename
fn = strrep(fn,'_ecog.mat','');

%% Define disease state
%Will output coh data into folders specific to PD, DYS, ET, or Epilepsy.
%The "Other" category for patients that do not fit well into
%standard categories (ie PD with tremor, PD with dystonic sx).
Dx = input('Enter patient diagnosis: 1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other : ');

%% Define ECOG contact over Motor Strip
% This allows for estimating M1-S1 coherence.
% M1_contact = input('Enter contact closest to M1 (format: [1 2 ...]): ');
% %-not necessary as M1_ch now imported from apmconv7 code

%% Define ecog contact over sensory cortex
% This allows for estimating M1-S1 coherence.
S1_ch = M1_ch -2;
if M1_ch==1
    S1_ch = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Function will terminate'],'OK');
    return
elseif M1_ch==2
    S1_ch = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Function will terminate'],'OK');
    return
end
%% Remontage Data
%such that contact 1 is referencing contact 2, 2 references 3, etc

num_contact_pair = length(ecog.contact_pair);
for i = 1:num_contact_pair  
    if i >= 5 
% the last ecog pair is already 5-6, thus does not require subtraction; last contact pair is actually the lfp 
            ecog.contact_pair(i).recog_signal = ecog.contact_pair(i).raw_ecog_signal; %recog=remontaged ecog
             % pt Solis ecog data, input ecog signal and reference contacts were
% accidentally flipped.  Correct for flipped ecog signal amplitude.
            %ecog.contact_pair(i).recog_signal = -ecog.contact_pair(i).raw_ecog_signal;
    else            
        %ecog.contact_pair(i).recog_signal = ecog.contact_pair(i).raw_ecog_signal - ...
                %ecog.contact_pair(i+1).raw_ecog_signal;
            % pt Solis ecog data, input ecog signal and reference contacts were
% accidentally flipped.  Correct for flipped ecog signal amplitude.
             ecog.contact_pair(i).recog_signal = ecog.contact_pair(i+1).raw_ecog_signal - ...
                ecog.contact_pair(i).raw_ecog_signal;
    end
end

%% Time-limited Parsing of rest/active epochs from each contact pair
% initialize structures that will contain all rest/active analysis
% This is specific to defined epoch lengths, rather than allowing analysis
% of the full length of each epoch, however long that may be
rest = struct('contact_pair',{});
active = struct('contact_pair',{});

num_epoch = length(ecog.active_time);

for i = 1:num_contact_pair

    % parse each rest/active epoch from each contact pair

    for j = 1:num_epoch
        start_rest = int32((ecog.rest_time(j) + OFFSET) * Fs); % time offset added to epoch times
%         end_rest = int32(start_rest + (Fs * EPOCH_LEN) - 1);
        end_rest = int32(ecog.active_time(j) * Fs - 1); %substract 1 because end of rest period is the datum just before the start of the active period that follows
        rest(1).contact_pair(i).epoch(j).recog_signal = ...
            ecog.contact_pair(i).recog_signal(start_rest:end_rest);

        start_active = int32((ecog.active_time(j)+ OFFSET) * Fs);
%         end_active = int32(start_active + Fs*EPOCH_LEN - 1);
        end_active = int32(ecog.rest_time(j+1) * Fs - 1); %substract 1 because end of rest period is the datum just before the start of the active period that follows
        active(1).contact_pair(i).epoch(j).recog_signal = ...
            ecog.contact_pair(i).recog_signal(start_active:end_active);
    end
end

%% Calculate coherence using mscohere
% perform coherence analysis using mscohere
% average all the epochs for rest/active

%ALC 09/07/2009: In some cases, will not have LFP data but still can
%calculate M1-S1 coherence. Code updated to do this. 

if num_contact_pair == 6
    for i = 1:num_epoch
        
        %       Comparing all remontaged contacts to the LFP
        for j = 1:num_contact_pair-1 % not including lfp
            
            % Calculate Ecog-LFP coherence
            inp1_r = rest.contact_pair(j).epoch(i).recog_signal;
            inp2_r = rest.contact_pair(6).epoch(i).recog_signal;
            
            inp1_a = active.contact_pair(j).epoch(i).recog_signal;
            inp2_a = active.contact_pair(6).epoch(i).recog_signal;
            
            %           [Cxy,F] = mscohere(x,y,window,noverlap,nfft,Fs)
            
            [Cxy,F] = mscohere(inp1_r,inp2_r,window,noverlap,nfft,Fs);
            all_epoch_rest(:,j,i) = Cxy; % store coh values in 3D matrix, coh x channels x epochs
            
            
            [Cxy,F] = mscohere(inp1_a,inp2_a,window,noverlap,nfft,Fs);
            all_epoch_active(:,j,i) = Cxy;
        end
    end
end
for i = 1:num_epoch
    % Calculate M1-S1 coherence
    inp1_r = rest.contact_pair(M1_ch).epoch(i).recog_signal;
    inp2_r = rest.contact_pair(S1_ch).epoch(i).recog_signal;
    
    inp1_a = active.contact_pair(M1_ch).epoch(i).recog_signal;
    inp2_a = active.contact_pair(S1_ch).epoch(i).recog_signal;
    
    %           [Cxy,F] = mscohere(x,y,window,noverlap,nfft,Fs)
    
    [Cxy,F] = mscohere(inp1_r,inp2_r,window,noverlap,nfft,Fs);
    all_epoch_M1S1rest(:,i) = Cxy; % store coh values in 3D matrix, coh x channels x epochs
    
    
    [Cxy,F] = mscohere(inp1_a,inp2_a,window,noverlap,nfft,Fs);
    all_epoch_M1S1active(:,i) = Cxy;
end

%Average all the epochs together. Resulting matrix will have 
%rows = window length, columns = comparisons (ecog contact vs LFP)

if num_contact_pair == 6
    rest.EcogLfp = mean(all_epoch_rest,3);
    active.EcogLfp = mean(all_epoch_active,3);
end
% M1S1.rest = mean(all_epoch_M1S1rest,2);
% M1S1.active = mean(all_epoch_M1S1active,2);
rest.M1S1 = mean(all_epoch_M1S1rest,2);
active.M1S1 = mean(all_epoch_M1S1active,2);
%% Calculate inverse hyperbolic tangent for transformed coherence 
% Y = atanh(X) returns the inverse hyperbolic tangent for each element of X.
if num_contact_pair == 6
    trans_coh.rest.EcogLfp = atanh(sqrt(rest.EcogLfp));
    trans_coh.active.EcogLfp = atanh(sqrt(active.EcogLfp));
end
% M1S1.trans_coh.rest = atanh(sqrt(M1S1.rest));
% M1S1.trans_coh.active = atanh(sqrt(M1S1.active));
trans_coh.rest.M1S1 = atanh(sqrt(rest.M1S1));
trans_coh.active.M1S1 = atanh(sqrt(active.M1S1));
trans_coh.freq = F;
trans_coh.M1 = M1_ch;

%% Save and write trans coh 
%Creates .mat output file of name (filename_transcoh.mat). transcoh.mat
%files can then be used to analyze coherence data across groups
if Dx==1
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\PD\'];
elseif Dx==2
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\DYS\'];
elseif Dx==3
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\ET\'];
elseif Dx==4
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\Epilepsy\'];
else
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\Other\'];
end
outputname = [outputdir fn,'_transcoh.mat'];
disp(['Writing ECOG/LFP  channels to:  ' outputname]);
save(outputname,'trans_coh');

%% Plotting Ecog-LFP Coherence
if num_contact_pair == 6
hf = figure('Name',(fn));

% coherence (0 to 1) vs F (0 to 150)
for i = 1:num_contact_pair-1
    h = subplot(4,5,i);
    plot(F, rest.EcogLfp(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
    hold on;
    plot(F, active.EcogLfp(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 30 0 0.5]) %axis of alpha and beta freq
    set (gca,'xtick',[0:20:150])
   if i==1
       ylabel('Coherence')
       xlabel('F')
   end
end
hl = legend('rest','active','Location','NorthWestOutside');
set(hl,'Position',[0.01016 0.8695 0.07109 0.0527]);

% coherence (0 to 0.5) vs F (0 to 50)
for i = 1:num_contact_pair-1 
    h = subplot(4,5,(i+5));
    plot(F, rest.EcogLfp(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
    hold on;
    plot(F, active.EcogLfp(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
    axis([0 50 0 0.5]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:0.1:0.5])
    if i==1
       ylabel('Coherence')
       xlabel('F')
    end
%     title(['Coherence of channel ',num2str(i),' and LFP']);
end

% Transformed coherence vs F
for i = 1:num_contact_pair-1
    h = subplot(4,5,(i+10));
    plot(F, trans_coh.rest.EcogLfp(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Transformed Coherence during Rest','color','b','LineWidth',1.5);
    hold on;
    plot(F, trans_coh.active.EcogLfp(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Transformed Coherence during Movement','color','r','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 30 0 0.5]) %axis of alpha and beta freq
    set (gca,'xtick',[0:20:150])
%     daspect([20 1 1])
    if i==1
       ylabel('Transformed Coherence')
       xlabel('F')
    end
%     set (gca,'xtick',[0:5:50])
%     set (gca,'ytick',[0:0.1:0.5])
%     title(['Transformed coherence of channel ',num2str(i),' and LFP']);
end

for i = 1:num_contact_pair-1 
    h = subplot(4,5,(i+15));
    plot(F, trans_coh.rest.EcogLfp(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
    hold on;
    plot(F, trans_coh.active.EcogLfp(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
    axis([0 50 0 1]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
    set (gca,'xtick',[0:10:50])
    set (gca,'ytick',[0:0.2:1])
%     daspect([10 1 1])
    if i==1
       ylabel('Transformed Coherence')
       xlabel('F')
    end
%     title(['Coherence of channel ',num2str(i),' and LFP']);
end

annotation(hf, 'textbox','String',fn,'HorizontalAlignment','left',...
    'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);

%% These annotations create the contact headers at top of graph
% Create textbox
annotation(hf,'textbox','String',{'C1 vs LFP'},'HorizontalAlignment','center','FontWeight','bold','FitHeightToText','off','LineStyle','none',...
    'Position',[0.1306 0.9222 0.1265 0.0485]);

% Create textbox
annotation(hf,'textbox','HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.2816 0.915 0.1265 0.05882]);

% Create textbox
annotation(hf,'textbox','String',{'C2 vs LFP'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.2867 0.9177 0.1265 0.05613]);

% Create textbox
annotation(hf,'textbox','String',{'C3 vs LFP'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.4486 0.927 0.1265 0.0485]);

% Create textbox
annotation(hf,'textbox','String',{'C4 vs LFP'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.6132 0.9177 0.1265 0.05613]);

% Create textbox
annotation(hf,'textbox','String',{'C5 vs LFP'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.7778 0.9301 0.1265 0.0485]);

%% Create rectangle around contact over motor strip

% M1_contact = input('Enter contact closest to M1 (format: [1 2 ...]): ');

if M1_ch == 1
    %if C1 over M1
    annotation(hf,'rectangle',[0.1109 0.0462 0.1609 0.9484]);
end

if M1_ch == 2
    %if C2 over M1
    annotation(hf,'rectangle',[0.2664 0.05027 0.1672 0.9484]);
end

if M1_ch == 3
    %if C3 over M1
    annotation(hf,'rectangle',[0.432 0.05027 0.1672 0.9484]);
end

if M1_ch == 4
    %if C4 over M1
    annotation(hf,'rectangle',[0.5914 0.04755 0.1672 0.9484]);
end

if M1_ch ==5
    %if C5 over M1
    annotation(hf,'rectangle',[0.757 0.0462 0.1672 0.9484]);
end

saveas(hf,[outputdir fn, '_transcoh'],'fig');
end
%% Plotting M1-S1 coherence

hf2 = figure('Name',(fn));
    plot(F, trans_coh.rest.M1S1, 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
    hold on;
    plot(F, trans_coh.active.M1S1, 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 30 0 0.5]) %axis of alpha and beta freq
    set (gca,'xtick',[0:20:150])
       ylabel('Transformed Coherence')
       xlabel('F')
       title('M1-S1 Coherence');

       annotation(hf2,'textbox','String',{fn},...
           'Position',[0.01135 0.9502 0.2789 0.03481],...
           'LineStyle','none', 'FitHeightToText','on');
        
clear