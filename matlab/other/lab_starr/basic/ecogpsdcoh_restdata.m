%this program graphs resting state PSD for ecog and lfp channels,
%and ecog-lfp coherence for each channel combination.  It should be run on
%resting activity where there are no movement-related epochs.  It requires
%file ecog_lfp_raw_data to be in the workspace. 

%first step: extract ecog bipolar recordings, and lfp recording, from
%ecog_lfp_raw_data

% use menu function to input yes or no to whether there is a meaningful lfp
% channel

%note: see mscript runspikeoscil to see robs code on converting timestamps
%to appropriate format (zeros and ones) to run mscohere

[filename pathname]=uigetfile('*.mat');
load([pathname filename])
filename=filename(1:end-13); % this takes off the ecog_lfp_mat ending of the filename;

ecog12=ecog_lfp_raw_data(:,1)-ecog_lfp_raw_data(:,2);
ecog23=ecog_lfp_raw_data(:,2)-ecog_lfp_raw_data(:,3);
ecog34=ecog_lfp_raw_data(:,3)-ecog_lfp_raw_data(:,4);
ecog45=ecog_lfp_raw_data(:,4)-ecog_lfp_raw_data(:,5);
ecog56=ecog_lfp_raw_data(:,5);
lfp=ecog_lfp_raw_data(:,6);

% %% Downsampling for Miller case
% ecog12 = downsample(ecog12,5);
% ecog23 = downsample(ecog23,5);
% ecog34 = downsample(ecog34,5);
% ecog45 = downsample(ecog45,5);
% ecog56 = downsample(ecog56,5);
% lfp = downsample(lfp,5);

%second step: generate power spectral density for each bipolar ecog
%recording and the lfp channel, and coherence spectrum for each ecog channel with lfp
%uses the Welch periodogram with 512 point fft for each

[psdecog12,F]=pwelch(ecog12,512,256,512,1000);
psdecog23=pwelch(ecog23,512,256,512,1000);
psdecog34=pwelch(ecog34,512,256,512,1000);
psdecog45=pwelch(ecog45,512,256,512,1000);
psdecog56=pwelch(ecog56,512,256,512,1000);
psdlfp=pwelch(lfp,512,256,512,1000);

cohecog12lfp=mscohere(ecog12,lfp,512,256,512,1000);
cohecog23lfp=mscohere(ecog23,lfp,512,256,512,1000);
cohecog34lfp=mscohere(ecog34,lfp,512,256,512,1000);
cohecog45lfp=mscohere(ecog45,lfp,512,256,512,1000);
cohecog56lfp=mscohere(ecog56,lfp,512,256,512,1000);

% Third step: plot the PSD's in first row with xscale 0-50 hz and y scale
% automatic, plot coherence in second row

psdall=[psdecog12 psdecog23 psdecog34 psdecog45 psdecog56 psdlfp];
cohall=[cohecog12lfp cohecog23lfp cohecog34lfp cohecog45lfp cohecog56lfp];
ymax=max(max(cohall));
ymax=ceil(10*ymax);
ymax=ymax/10;
for i=1:11
    if i<=6
        subplot(2,6,i)
        plot(F,psdall(:,i))
        xlim([0 50]);
        if i==1;
            str=[filename sprintf('\n') 'psdecog' num2str(i) num2str(i+1)]; % allows title to have file name
        elseif i==6
            str='lfp';
        else
            str=['psdecog' num2str(i) num2str(i+1)];
        end
        title(str);
    else
        subplot(2,6,i)
        plot(F,cohall(:,i-6)); %note that cohall starts with column 1 not column 7 thus i-6 is used
        xlim([0 50]);
        ylim([0 ymax]);
        str=['cohecog' num2str(i-6) num2str(i-5) 'lfp'];
        title(str);
    end
end

        
        


