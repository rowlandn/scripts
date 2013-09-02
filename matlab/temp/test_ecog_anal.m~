%% practice ECOG data analysis

% on next run, describe inputs and outputs of each script

% addpaths
addpath ~/Dropbox/cluster_files/data_raw/recordings/ecog/low_res/bruner
addpath ~/Dropbox/cluster_files/libraries/matlab/scripts/other/lab_starr/mat_scripts/alphaomega;
addpath ~/Dropbox/cluster_files/libraries/matlab/scripts/other/lab_starr/mat_scripts/both;
addpath ~/Dropbox/cluster_files/libraries/matlab/scripts/other/lab_starr/mat_scripts/mat_for_high_res_grid;

% examine SSEP
cd ~/Dropbox/cluster_files/data_raw/recordings/ecog/low_res/bruner

load BrunL01LecogSSEP

plot(CECOG_1,'b'); hold on
plot(CECOG_2,'r')
plot(CECOG_3,'k')
plot(CECOG_4,'c')
plot(CECOG_5,'m') % next time, subtract the contacts and see if you can see this is raw data

aomatSSEP % determine M1 contact

% converts .json file into .mat file with descriptor strings and 
% timestamp values for each phase
% input = .json file
% output = .mat file with 2 variables: description string and timestamp
% numeric values

import_from_ipad

% converts mat file from AlphaOmega Mapfile converter program to a format
% that can be run on existing code for ecogPSD analysis and new ipad
% analyses (essentially aligns and scales all channels then places into
% common structure)
%
% input: raw .mat file of all channels during iPad task + M1 channel
% output: .mat file of same name and '_ecog' appended

clear

aomatconv_ipad('~/Dropbox/cluster_files/data_raw/recordings/ecog/low_res/bruner/BrunL15LecogLlfp_postlead_iPad',3);

plot(aux.chan(1).raw)
plot(aux.chan(2).raw)
plot(aux.chan(3).raw)
plot(emg.chan(1).raw)
plot(emg.chan(2).raw)
plot(emg.chan(3).raw)
plot(ecog.contact_pair(1).raw_ecog_signal)
plot(ecog.contact_pair(2).raw_ecog_signal)
plot(ecog.contact_pair(3).raw_ecog_signal)

%load BrunL15LecogLlfp_postlead_iPad_ecog



% Determine epochs visually(?)

Determine_events

%% plot psd's


timePSD_ipad_PREP
timePSD_ipad_MVT


