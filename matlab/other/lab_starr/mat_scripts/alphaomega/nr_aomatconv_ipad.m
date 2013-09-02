function aomatconv_ipad(filename,M1_ch)
% converts mat file from AlphaOmega Mapfile converter program to a format
% that can be run on existing code for ecogPSD analysis and new ipad
% analyses

% Created by SAS (12/15/2008)pl
% Edited by SAS (5/29/09) to process EMG/accel data.
% 
% Input:
%       1) "fname.mat" = output file form Mapfileconverter
% Output:
%       1) "fname_ecog.mat" = ecog/LFP data


%% load data

% [filename pathname]=uigetfile('*.mat','Select .mat file');
% cd(pathname);

load(filename);
filename= strrep(filename,'.mat','');
assignin('base','filename',filename)
%assignin('base','filename',filename)
% make sure sampling frequencies are equal for ecog/LFP/task recordings, then
% use CECOG_1_KHz varible as the sampling rate
% if isequal(CECOG_1_KHz,CECOG_2_KHz,CECOG_3_KHz,CECOG_4_KHz,CECOG_5_KHz, CLFP_KHz,CAn_In__1_KHz,CAn_In__2_KHz,CAn_In__3_KHz,CEMG_1_KHz,CEMG_2_KHz,CEMG_3_KHz)
if exist('CECOG1_KHz')
    Fs = CECOG1_KHz;
    Fs = Fs*1000;
else
    Fs = CECOG_1_KHz;
    Fs = Fs*1000; 
end
% multiply by 1000 for KHz->Hz conversion
% else
%     error('Sampling frequencies are different for ecog/LFP data');
% end

%% M1 contact
% M1_ch = input('Enter channel # for M1 contact (usually 4 or 5): '); % added 9/9/09 to be carried over for later programs

%% force ecog/lfp data to be the same lengths for montage-ing

minlength=min([length(CECOG_1) length(CECOG_2) length(CECOG_3)...
    length(CECOG_4) length(CECOG_5) length(CLFP)]);
CECOG_1  = CECOG_1(1:minlength);
CECOG_2  = CECOG_2(1:minlength);
CECOG_3  = CECOG_3(1:minlength);
CECOG_4  = CECOG_4(1:minlength);
CECOG_5  = CECOG_5(1:minlength);
CLFP   = CLFP(1:minlength);

%% ecog/LFP voltage calibration
C1 = 5/32768; % constant for calibration
% note: value 32768 is equal to 5V in the Alpha Omega system
C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Gecog = 7000;   % ecog channel gain
Glfp = 25000;   % LFP channel gain

CAn_In__1 = CAn_In__1*C1;
CAn_In__2 = CAn_In__2*C1;
CAn_In__3 = CAn_In__3*C1;

% CEMG_1 = CEMG_1*C1;
% CEMG_2 = CEMG_2*C1;
% CEMG_3 = CEMG_3*C1;
% 
CECOG_1  = CECOG_1*C1*C2/Gecog;
CECOG_2  = CECOG_2*C1*C2/Gecog;
CECOG_3  = CECOG_3*C1*C2/Gecog;
CECOG_4  = CECOG_4*C1*C2/Gecog;
CECOG_5  = CECOG_5*C1*C2/Gecog;
CLFP   = CLFP*C1*C2/Glfp;

%% store processed ecog/LFP data into ecog structure array

% frq1 = menu('Select one','resample1K','resample2K','no');

d = round(Fs/1000);
if int32(Fs) == 1502
    frq1 = 2;
    frq2 = 3;
    d=Fs/1000;
else
    frq1=1;
    frq2=d;
    d=Fs/1000;
    
end

% ecog = struct('contact_pair',{},'rest_time',{},'active_time',{});
% evt = menu('Select one','LFP','No LFP');

ecog(1).contact_pair(1).raw_ecog_signal=resample(CECOG_1,frq1,frq2);
ecog(1).contact_pair(2).raw_ecog_signal=resample(CECOG_2,frq1,frq2);
ecog(1).contact_pair(3).raw_ecog_signal=resample(CECOG_3,frq1,frq2);
ecog(1).contact_pair(4).raw_ecog_signal=resample(CECOG_4,frq1,frq2);
ecog(1).contact_pair(5).raw_ecog_signal=resample(CECOG_5,frq1,frq2);

if ~isempty(strfind(filename,'lfp')) || ~isempty(strfind(filename,'LFP'))% evt == 1
    ecog(1).contact_pair(6).raw_ecog_signal=resample(CLFP,frq1,frq2);
end
Fs = 1000;

%% store processed Aux chan data into aux structure array

aux = struct('chan',{});
aux(1).chan(1).raw=resample(CAn_In__1,frq1,frq2);
aux(1).chan(2).raw=resample(CAn_In__2,frq1,frq2);
aux(1).chan(3).raw=resample(CAn_In__3,frq1,frq2);

%% store processed EMG chan data into aux structure array

emg = struct('chan',{});
emg(1).chan(1).raw=resample(CEMG_1,frq1,frq2);
emg(1).chan(2).raw=resample(CEMG_2,frq1,frq2);
emg(1).chan(3).raw=resample(CEMG_3,frq1,frq2);

%% remove zeros add when map data are converted in mat file by 'mapconverter'


%% Do the remontage
for i = 1: length(ecog.contact_pair)
    if i<5
        ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal - ecog.contact_pair(i+1).raw_ecog_signal;
    else
        ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal;
    end
end

%% grab rest/active timestamps from task button

% ask user for aux channel with task button voltage
% evt = menu('Select one','rest/active data','rest only');

if ~isempty(strfind(filename,'ipad')) %evt == 1
    
    % find active and rest periods using the raw signal
    MARGIN = 0.5; % is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
    thresh = max(CAn_In__3); 
    ACTIVE = [thresh-MARGIN thresh+MARGIN]; %Determine the threshold
    in1 = inrange(CAn_In__3,ACTIVE);
    inds = find(in1);
    [pos,n] = evFindGroups(inds,1,1000*d); %find active period of 1s minimum
    active_time = inds(pos(1,:))/d;
    rest_time = inds(pos(2,:))/d;
    
    %assignin('base','thresh',thresh)
    %%assignin('base','active',active)
    %assignin('base','in1',in1)
    %%assignin('base','inds',inds)
    %assignin('base','active_time',active_time)
    %assignin('base','rest_time',rest_time)
        
    % save rest/active timestamps in ecog structure array
    ecog(1).rest_time=int32(rest_time);
    ecog(1).active_time=int32(active_time);
    
    % find the begining of each trial using the raw signal
    CAn_In__1=CAn_In__1-mean(CAn_In__1);
%     thresh = max(CAn_In__1);
%     MARGIN = thresh/5; % 0.3 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
%     ACTIVE = [thresh-MARGIN thresh+MARGIN];
%     in1 = inrange(CAn_In__1,ACTIVE);
%     inds = find(in1);
%     thresh = max(CAn_In__1)-3*std(CAn_In__1);
    
    close ;plot(CAn_In__1)
    hold on
    thresh =  0.02 %input('thresh');
    
    %assignin('base','CAn_In__1',CAn_In__1)
    %assignin('base','thresh',thresh)
    
     
    inds = find(abs(CAn_In__1)>thresh);
    
    assignin('base','inds',inds)
    
    [pos,n] = evFindGroups(inds,5*d,2);
    ipad_ON_time = inds(pos(1,1:end-1))/d;
    ipad_OFF_time = inds(pos(1,2:end))/d;
    % delete false detection
    diff_ipad_ON=  ipad_ON_time(2:end)-ipad_ON_time(1:end-1);
    ok_ON = find(diff_ipad_ON>5000);
    
    %%assignin('base','evFindGroups',evFindGroups)
    %assignin('base','d',d)
    %assignin('base','thresh',thresh)
    %%assignin('base','active',active)
    %assignin('base','in1',in1)
    %%assignin('base','inds',inds)
    %assignin('base','active_time',active_time)
    %assignin('base','rest_time',rest_time)
    
    
    
    
    
    %assignin('base','ipad_ON_time',ipad_ON_time)
    %assignin('base','diff_ipad_ON',diff_ipad_ON)
    %assignin('base','ok_ON',ok_ON)
    
    diff_ipad_OFF=  ipad_OFF_time(2:end)-ipad_OFF_time(1:end-1);
    ok_OFF = find(diff_ipad_OFF>5000);
    % check the timing
    plot(ipad_ON_time(ok_ON)*d, thresh,'*r')
    plot(ipad_OFF_time(ok_OFF)*d, thresh,'*k')
     
    ecog(1).ipad_ON_time=int32(ipad_ON_time(ok_ON));
    ecog(1).ipad_OFF_time= int32(ipad_OFF_time(ok_OFF));
    
    %assignin('base','ok_ON',ok_ON)
    %assignin('base','ipad_ON_time',ipad_ON_time)

else
    ecog(1).rest_time=[];
    ecog(1).active_time=[];
    ecog(1).ipad_ON_time=[];
    ecog(1).ipad_OFF_time=[];
    
end

%% save data
% use filename to create output name

name=[filename(1:11),'anl_rec_',filename(12:end)];

assignin('base','name',name)

cd(['~/Dropbox/cluster_files/data/',filename(1:10),'/anl'])
save(name,'name','ecog','M1_ch','aux','emg','Fs');

%% Add IPAD file info
[ipad_data_filename, ipad_data_dir] = uigetfile('*.mat','Choose ipad data file'); % choose ipad data file
load(ipad_data_filename)

%assignin('base','ipad_data_dir',ipad_data_dir)
% 
% 
% name_ipad=[filename '_ecog'];
% 
% load(name_ipad);
%cd([ipad_data_dir(1:end-5),'anl'])


cd(['~/Dropbox/cluster_files/data/',filename(1:10),'/anl'])
save(name,'description','timestamp','-append');