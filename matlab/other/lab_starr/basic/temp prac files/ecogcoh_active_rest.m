%% Creating separate bins for active and rest data
% 8/28/08 Got stuck on num_contact_pair = length(ecog.contact_pair) with
% error message: ecog not a recognized function or variable

%%
% import filename and pathname of mat file containing ecog data created by
% apmconv7_coh
[fn pn] = uigetfile('*.mat','Select .mat containing ecog_lfp data');
cd(pn);
% load ecog_rest and ecog_active variables and store both into one
% structure called ecog_lfp
load([pn fn]);
% remove '_ecog_lfp.mat' ending from filename
fn = strrep(fn,'_ecog_lfp.mat','');

hf = figure;
num_contact_pair = length(ecog_lfp.contact_pair);
num_epoch = length(ecog_lfp.rest_time);

% initialize structures that will contain all rest/active analysis
rest = struct('contact_pair',{});
active = struct('contact_pair',{});

           
for i = 1:num_contact_pair
    
    % parse each rest/active epoch from each contact pair
    
    for j = 1:num_epoch
        start_rest = (ecog.rest_time(j) + OFFSET) * Fs; % time offset added to epoch times
        end_rest = start_rest + (Fs * EPOCH_LEN) - 1; 
        rest(1).contact_pair(i).epoch(j).raw_ecog_signal = ...
            ecog.contact_pair(i).raw_ecog_signal(start_rest:end_rest);
        start_active = (ecog.active_time(j)+ OFFSET) * Fs;
        end_active = start_active + Fs*EPOCH_LEN - 1;
        active(1).contact_pair(i).epoch(j).raw_ecog_signal = ...
            ecog.contact_pair(i).raw_ecog_signal(start_active:end_active);
    end
    
end