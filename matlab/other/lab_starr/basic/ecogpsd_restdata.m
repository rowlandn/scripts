%this program graphs resting state PSD for ecog channels
%resting activity where there are no movement-related epochs.  It requires
%a matfile containing ecog_lfp_raw data and loads this into the workspace
 

%first step: extract ecog bipolar recordings from
%ecog_lfp_raw_data


[filename pathname]=uigetfile('*.mat');
load([pathname filename])
filename=filename(1:end-13); % this takes off the ecog_lfp_mat ending of the filename;

ecog12=ecog_lfp_raw_data(:,1); %-ecog_lfp_raw_data(:,2);
ecog23=ecog_lfp_raw_data(:,2); %-ecog_lfp_raw_data(:,3);
ecog34=ecog_lfp_raw_data(:,3); %-ecog_lfp_raw_data(:,4);
ecog45=ecog_lfp_raw_data(:,4); %-ecog_lfp_raw_data(:,5);
ecog56=ecog_lfp_raw_data(:,5);


%second step: generate power spectral density for each bipolar ecog
%recording uses the Welch periodogram with 512 point fft for each

[psdecog12,F]=pwelch(ecog12,512,256,512,1000);
psdecog23=pwelch(ecog23,512,256,512,1000);
psdecog34=pwelch(ecog34,512,256,512,1000);
psdecog45=pwelch(ecog45,512,256,512,1000);
psdecog56=pwelch(ecog56,512,256,512,1000);




% Third step: plot the PSD's in first row with xscale 0-50 hz and y scale
% automatic, plot coherence in second row

psdall=[psdecog12 psdecog23 psdecog34 psdecog45 psdecog56];
for i=1:5
        subplot(1,5,i)
        plot(F,psdall(:,i))
        xlim([0 50]);
        if i==1;
            str=[filename sprintf('\n') 'psdecog' num2str(i) num2str(i+1)]; % allows title to have file name
        else
            str=['psdecog' num2str(i) num2str(i+1)];
        end
        title(str);
end

        
        


