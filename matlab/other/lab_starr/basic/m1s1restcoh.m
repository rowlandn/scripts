%% Define Constants

OFFSET = 0.5; % epoch time offset in seconds, NOTE we may change this to zero
Fs=1000;            % samples per unit time, in this case digitization at 1000 hz rate
nfft = 512; %FFT length that determines the frequencies at which the coherence is estimated. For real x and y,the length of Cxy is (nfft/2+1) if nfft is even or (nfft+1)/2 if nfft is odd. 
window = 512; %epoch data will be divided into chunks of size equal to window for coherence calculation
noverlap = 256; % how much overlap there is between windows
% nfft and window must be equal and noverlap is a fraction of those (50%,
% 75%, etc.


%% Import/Prepare file
% import filename and pathname of mat file containing ecog data created by
% apmconv7_coh
[fn pn] = uigetfile('*.mat','Select .mat containing ecog_lfp data');
cd(pn);
load([pn fn]);
% remove '_ecog_lfp.mat' ending from filename
fn = strrep(fn,'_ecog_lfp.mat','');

%% Define disease state
%Will output coh data into folders specific to PD or DYS
%Code not currently flexible enough to handle alterate diagnoses(ALC 11/25/08)
Dx = input('Enter patient diagnosis: 1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other : ');

% switch Dx  --attempt to make error message if wrong number entered
%     case
%         disp([num2str(input) 'is not a valid diagnosis']);
% end

%% Define ECOG contact over Motor Strip
M1_contact = input('Enter contact closest to M1 (format: [1 2 ...]): ');

%% Define ecog contact over sensory cortex
S1_contact = M1_contact -2;
if M1_contact==1
    s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Function will terminate'],'OK');
    return
elseif M1_contact==2
    s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Function will terminate'],'OK');
    return
end
%% Remontage Data
%such that contact 1 is referencing contact 2, 2 references 3, etc

% Recog = zeros(size(ecog_lfp_raw_data(:,end-2);
ecog1v2 = ecog_lfp_raw_data(:,1) - ecog_lfp_raw_data(:,2);
ecog2v3 = ecog_lfp_raw_data(:,2) - ecog_lfp_raw_data(:,3);
ecog3v4 = ecog_lfp_raw_data(:,3) - ecog_lfp_raw_data(:,4);
ecog4v5 = ecog_lfp_raw_data(:,4) - ecog_lfp_raw_data(:,5);
ecog5v6 = ecog_lfp_raw_data(:,5);
LFP = ecog_lfp_raw_data(:,end-2);
% %% Downsampling for Miller case (accidentally recorded at 5K)
% ecog1v2 = downsample(ecog1v2,5);
% ecog2v3 = downsample(ecog2v3,5);
% ecog3v4 = downsample(ecog3v4,5);
% ecog4v5 = downsample(ecog4v5,5);
% ecog5v6 = downsample(ecog5v6,5);
% LFP = downsample(LFP,5);

recog = [ecog1v2 ecog2v3 ecog3v4 ecog4v5 ecog5v6 LFP];


%% Now add coherence

inp1_r = recog(:,int8(M1_contact));
inp2_r = recog(:,int8(S1_contact));
            
[Cxy,F] = mscohere(inp1_r,inp2_r,window,noverlap,nfft,Fs);
m1s1coh = Cxy;
 
%% Calculate inverse hyperbolic tangent for transformed coherence 
% Y = atanh(X) returns the inverse hyperbolic tangent for each element of X.

trans_coh.m1s1.trurest = atanh(sqrt(m1s1coh));
trans_coh.freq = F;

%% Store trans coh data for contact over M1
%This will create separate fields for storing coh data from directly over
%M1. This will be handy in later code for averaging coherence between motor
%cortex and STN across subjects
trans_coh.M1 = M1_contact

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
outputname = [outputdir fn,'_m1s1coh.mat'];
disp(['Writing ECOG/LFP  channels to:  ' outputname]);
save(outputname,'trans_coh');

%% Plotting

hf = figure('Name',(fn));
    plot(F, trans_coh.m1s1.trurest, 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 30 0 0.5]) %axis of alpha and beta freq
    set (gca,'xtick',[0:20:150])
    ylabel('Coherence')
    xlabel('F')
    title('M1-S1 Coherence in Rest Only condition');
    annotation(hf,'textbox','String',{fn},...
           'Position',[0.01135 0.9502 0.2789 0.03481],...
           'LineStyle','none', 'FitHeightToText','on');

