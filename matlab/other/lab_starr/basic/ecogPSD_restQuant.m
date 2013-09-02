%this program graphs resting state PSD for ecog channels
%resting activity where there are no movement-related epochs.  It requires
%a matfile containing ecog_lfp_raw data and loads this into the workspace
 
% 4/9/09: updated to also include quantPSD and exportdata:ALC

%% Variables
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   
%% Get file
[filename pathname]=uigetfile('*.mat');
load([pathname filename])
filename = strrep(filename,'_ecog_lfp.mat',''); % this takes off the ecog_lfp_mat ending of the filename;
%% Define disease state
%Will output data into folders specific to disease state
Dx = input('Enter patient diagnosis: 1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other : ');

%% first step: extract ecog bipolar recordings from
%ecog_lfp_raw_data
ecog12=ecog_lfp_raw_data(:,1); %-ecog_lfp_raw_data(:,2);
ecog23=ecog_lfp_raw_data(:,2); %-ecog_lfp_raw_data(:,3);
ecog34=ecog_lfp_raw_data(:,3); %-ecog_lfp_raw_data(:,4);
ecog45=ecog_lfp_raw_data(:,4); %-ecog_lfp_raw_data(:,5);
ecog56=ecog_lfp_raw_data(:,5);
lfp=ecog_lfp_raw_data(:,6);

%second step: generate power spectral density for each bipolar ecog
%recording uses the Welch periodogram with 512 point fft for each

[psdecog12,F]=pwelch(ecog12,512,256,512,1000);
psdecog23=pwelch(ecog23,512,256,512,1000);
psdecog34=pwelch(ecog34,512,256,512,1000);
psdecog45=pwelch(ecog45,512,256,512,1000);
psdecog56=pwelch(ecog56,512,256,512,1000);
psdlfp=pwelch(lfp,512,256,512,1000);



% Third step: plot the PSD's in first row with xscale 0-50 hz and y scale
% automatic, plot coherence in second row

psdall=[psdecog12 psdecog23 psdecog34 psdecog45 psdecog56 psdlfp];
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

%% quantitative analysis of PSD
% now that all the data crunching is done and the results are stored in
% a 'psdall' matrix, run sub-function quantPSDrest for
% quantitative analysis of PSD data

[allfreq subfreq order] = quantPSDrest(psdall,FREQ_QPSD,F,filename);
        
%% Export data to excel
[status, message] = exportRESTdata(allfreq, subfreq, filename);
cd(pathname); %exportdata changes the path in order to write the excel data. This brings back current path.

%% Save data
if Dx==1
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\PD\rest'];
elseif Dx==2
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\DYS\rest'];
elseif Dx==3
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\ET\rest'];
elseif Dx==4
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\Epilepsy\'];
else
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\Other\rest'];
end
outputname = [outputdir filename,'_ecogPSD.mat'];
disp(['Writing ecog PSD data to:  ' outputname]);
% save(outputname,'ecogPSD');
save(outputname,'order','psdall', 'F', 'allfreq', 'subfreq');


