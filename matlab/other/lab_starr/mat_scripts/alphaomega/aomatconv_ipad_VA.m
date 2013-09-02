function aomatconv_ipad_VA(filename,M1_ch)
% converts mat file from AlphaOmega Mapfile converter program to a format
% that can be run on existing code for ecogPSD analysis and new ipad
% analyses

% Created by SAS (12/15/2008)
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


ipad_signal=CAN_IN2;
trig_signal=CAN_IN1;
% make sure sampling frequencies are equal for ecog/LFP/task recordings, then
% use CECOG1_KHz varible as the sampling rate
% if isequal(CECOG1_KHz,CECOG2_KHz,CECOG3_KHz,CECOG4_KHz,CECOG5_KHz, CLFP1_KHz,CAN_IN1_KHz,CAN_IN2_KHz,CAN_IN3_KHz,CEMG1_KHz,CEMG2_KHz,CEMG3_KHz)
    Fs = CECOG1_KHz;
    Fs = Fs*1000;                           % multiply by 1000 for KHz->Hz conversion
% else
%     error('Sampling frequencies are different for ecog/LFP data');
% end

%% M1 contact
% M1_ch = input('Enter channel # for M1 contact (usually 4 or 5): '); % added 9/9/09 to be carried over for later programs

%% force ecog/lfp data to be the same lengths for montage-ing
% 
% minlength=min([length(CECOG1) length(CECOG2) length(CECOG3)...
%     length(CECOG4) length(CECOG5) length(CLFP1)]);
% CECOG1  = CECOG1(1:minlength);
% CECOG2  = CECOG2(1:minlength);
% CECOG3  = CECOG3(1:minlength);
% CECOG4  = CECOG4(1:minlength);
% CECOG5  = CECOG5(1:minlength);
% CLFP1   = CLFP1(1:minlength);

%% ecog/LFP voltage calibration
C1 = 5/32768; % constant for calibration
% note: value 32768 is equal to 5V in the Alpha Omega system
C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Gecog = 7000;   % ecog channel gain
Glfp = 25000;   % LFP channel gain

CAN_IN1 = CAN_IN1*C1;
CAN_IN2 = CAN_IN2*C1;
CAN_IN3 = CAN_IN3*C1;

% CEMG1 = CEMG1*C1;
% CEMG2 = CEMG2*C1;
% CEMG3 = CEMG3*C1;
% 
CECOG1  = CECOG1*C1*C2/Gecog;
CECOG2  = CECOG2*C1*C2/Gecog;
CECOG3  = CECOG3*C1*C2/Gecog;
CECOG4  = CECOG4*C1*C2/Gecog;
CECOG5  = CECOG5*C1*C2/Gecog;
CLFP1   = CLFP1*C1*C2/Glfp;

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

ecog(1).contact_pair(1).raw_ecog_signal=resample(CECOG1,frq1,frq2);
ecog(1).contact_pair(2).raw_ecog_signal=resample(CECOG2,frq1,frq2);
ecog(1).contact_pair(3).raw_ecog_signal=resample(CECOG3,frq1,frq2);
ecog(1).contact_pair(4).raw_ecog_signal=resample(CECOG4,frq1,frq2);
ecog(1).contact_pair(5).raw_ecog_signal=resample(CECOG5,frq1,frq2);

if ~isempty(strfind(filename,'lfp')) || ~isempty(strfind(filename,'LFP'))% evt == 1
    ecog(1).contact_pair(6).raw_ecog_signal=resample(CLFP1,frq1,frq2);
end
Fs = 1000;

%% store processed Aux chan data into aux structure array

aux = struct('chan',{});
aux(1).chan(1).raw=resample(CAN_IN1,frq1,frq2);
aux(1).chan(2).raw=resample(CAN_IN2,frq1,frq2);
aux(1).chan(3).raw=resample(CAN_IN3,frq1,frq2);

%% store processed EMG chan data into aux structure array

emg = struct('chan',{});
emg(1).chan(1).raw=resample(CEMG1,frq1,frq2);
emg(1).chan(2).raw=resample(CEMG2,frq1,frq2);
emg(1).chan(3).raw=resample(CEMG3,frq1,frq2);


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
    thresh = max(trig_signal); 
    ACTIVE = [thresh-MARGIN thresh+MARGIN]; %Determine the threshold
    in1 = inrange(trig_signal,ACTIVE);
    inds = find(in1);
    [pos,n] = evFindGroups(inds,1,1000*d); %find active period of 1s minimum
    active_time = inds(pos(1,:))/d;
    rest_time = inds(pos(2,:))/d;
        
    % save rest/active timestamps in ecog structure array
    ecog(1).rest_time=int32(rest_time);
    ecog(1).active_time=int32(active_time);
    
    % find the begining of each trial using the raw signal
%     thresh =  max(CAN_IN1);
%     MARGIN = thresh/5; % 0.3 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
%     ACTIVE = [thresh-MARGIN thresh+MARGIN];
%     in1 = CAN_IN1>0.015;%inrange(CAN_IN1,ACTIVE);
%     inds = find(in1);
%     [pos,n] = evFindGroups(inds,100*d,2*d);%evFindGroups(inds,5*d,2*d);
    plot(ipad_signal)
    thresh =  input('thresh');
    close 
    y=abs(hilbert(ipad_signal));
    inds = find(abs(ipad_signal)>thresh);
    [pos,n] = evFindGroups(inds,5*d,2);
    ipad_ON_time = inds(pos(1,1:end-1))/d;
    ipad_OFF_time = inds(pos(1,2:end))/d;
    % delete false detection
    diff_ipad_ON=  ipad_ON_time(2:end)-ipad_ON_time(1:end-1);
    ok_ON = find(diff_ipad_ON>10000);
    diff_ipad_OFF=  ipad_OFF_time(2:end)-ipad_OFF_time(1:end-1);
    ok_OFF = find(diff_ipad_OFF>10000);
    
    ecog(1).ipad_ON_time=int32(ipad_ON_time(ok_ON));
    ecog(1).ipad_OFF_time= int32(ipad_OFF_time(ok_OFF));

else
    ecog(1).rest_time=[];
    ecog(1).active_time=[];
    ecog(1).ipad_ON_time=[];
    ecog(1).ipad_OFF_time=[];
    
end

%% save data
% use filename to create output name

name=[filename(1:end-4),'_ecog'];
save(name,'name','ecog','M1_ch','aux','emg','Fs');

%% Add IPAD file info
name_ipad=strrep(filename,'.mat','_ipad.mat');

try 
    load(name_ipad);
    save(name,'description','timestamp','-append');
end

