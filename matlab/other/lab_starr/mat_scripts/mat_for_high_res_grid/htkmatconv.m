
function htkmatconv(filename,M1_ch1,M1_ch2)

ecog = struct('chan',{});
aux = struct('chan',{});
aux(1).chan(1).raw=[];
aux(1).chan(2).raw=[];
aux(1).chan(3).raw=[];
emg = struct('chan',{});
X = struct('chan',{});
signal = struct('chan',{});

gain=1e6;
Fs=1000; % sampling freq after downsampling

files = dir('*.htk');


%% load the data and store processed Aux chan data into structure array
for abc =1 : length(files)
    if ~isempty(strfind(files(abc).name,'htk'))
        % load the data and convert in mat.file
        file_name=files(abc).name;
        [d,fs,dt,tc,t]=readhtk(file_name);
       
        %store the downsampled data in structure arrays
        if ~isempty(strfind(file_name,'ipad.htk'))
            aux(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'accel.htk'))
            aux(1).chan(2).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'trig.htk'))
            aux(1).chan(3).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'emg.htk'))
            emg(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'emg2.htk'))
            emg(1).chan(2).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'signal.htk'))
            signal(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        else            
            ecog(1).contact_pair(tc).raw_ecog_signal=resample(d,2^10,5^5)*gain;
        end
        abc
        length(ecog(1).contact_pair)
    end
end

C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Glfp = 25000;   % LFP channel gain
if ~isempty(signal)
    ecog(1).contact_pair(29).raw_ecog_signal=signal(1).chan(1).raw*C2/Glfp;
end

%% notch filter around 60Hz, 120Hz and 180Hz
% butterworth notch filter - model order, [low/(Fs/2) high/(Fs/2)]
[n1_b, n1_a]=butter(3,2*[57 63]/Fs,'stop'); %60hz
[n2_b, n2_a]=butter(3,2*[117 123]/Fs,'stop'); %120hz
[n3_b, n3_a]=butter(3,2*[177 183]/Fs,'stop'); %180hz
[n4_b, n4_a]=butter(3,2*[237 243]/Fs,'stop'); %240hz
[n5_b, n5_a]=butter(3,2*[297 303]/Fs,'stop'); %300hz
[n6_b, n6_a]=butter(3,2*[357 363]/Fs,'stop'); %360hz
%     [n7_b, n7_a]=butter(3,2*[88 94]/Fs,'stop'); %180hz
for k=1:length(ecog(1).contact_pair)
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n1_b, n1_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n2_b, n2_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 120
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n3_b, n3_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n4_b, n4_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n5_b, n5_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 120
    ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n6_b, n6_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
    %          ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n7_b, n7_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
end
%% remove DC offset
for k=1:length(ecog(1).contact_pair)
    ecog(1).contact_pair(k).raw_ecog_signal=ecog(1).contact_pair(k).raw_ecog_signal-mean(ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
end
for k=1:length(aux.chan)
    aux.chan(k).raw=aux.chan(k).raw-mean(aux.chan(k).raw);
end
% %% contact to remove
% 
% unused=input('unused contacts');
% used = setdiff(1:length(ecog.contact_pair),unused);
% for i = 1:length(used)
%     j=used(i);
%     X(1).contact_pair(i).raw_ecog_signal=ecog.contact_pair(j).raw_ecog_signal;
% end
% ecog=X;

%% detection of bad electrodes
% figure;hold on
% for i = 1: 14
%     subplot(3,5,i)
%     plot(ecog.contact_pair(i).raw_ecog_signal)
%     title(num2str(i))
% %     ylim([-0.005 0.005])
% end
% figure;hold on
% for i = 15: 28
%     subplot(3,5,i-14)
%     plot(ecog.contact_pair(i).raw_ecog_signal)
%     title(num2str(i))
% %     ylim([-25 25])
% end
% bad=input('bad contacts');
% 
% close all
%% common reference
if length(ecog.contact_pair) <=32
     length_CAR = 28;
elseif length(ecog.contact_pair) >=64
    length_CAR = 64;
end
data = nan*ones(length_CAR,length(ecog.contact_pair(1).raw_ecog_signal));
for i = 1: length_CAR
    data(i,:) = ecog.contact_pair(i).raw_ecog_signal';
end
car=nanmean(data)/length_CAR;
for i = 1: length_CAR
    ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal-car;
%     ecog.contact_pair(i).remontaged_ecog_signal=ecog.contact_pair(i).raw_ecog_signal-ecog.contact_pair(i+1).raw_ecog_signal;
end
% ecog.contact_pair(28).remontaged_ecog_signal=ecog.contact_pair(28).raw_ecog_signal;
if size(ecog.contact_pair,2) ==29
    ecog.contact_pair(29).remontaged_ecog_signal=ecog.contact_pair(29).raw_ecog_signal;
elseif size(ecog.contact_pair,2) >28 && size(ecog.contact_pair,2) <=32
    for ii = 29:size(ecog.contact_pair,2)-1
        ecog.contact_pair(ii).remontaged_ecog_signal=ecog.contact_pair(ii).raw_ecog_signal-ecog.contact_pair(ii+1).raw_ecog_signal;
    end
end
%% load ipad data

if ~isempty(strfind(filename,'_ipad')) %evt == 1
    
%     [d,fs,dt,tc,t]=readhtk('trig.htk');
%     chan_trig=d;
%     % find active and rest periods using the raw signal
%     MARGIN = 0.005; % is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
%     
% %     thresh = max(chan_trig);
% %     ACTIVE = [thresh-MARGIN thresh+MARGIN]; %Determine the threshold
% %     in1 = inrange(chan_trig,ACTIVE);
%     thresh =  4;%mean(chan_trig(find(chan_trig~=0)))
%     in1=chan_trig>thresh;
%     inds = find(in1);
%     [pos,n] = evFindGroups(inds,1,1000); %find active period of 1s minimum
%     active_time = inds(pos(1,:))/(fs/1000);
%     rest_time = inds(pos(2,:))/(fs/1000);
%     
%     % save rest/active timestamps in ecog structure array
%     ecog(1).rest_time=int32(rest_time);
%     ecog(1).active_time=int32(active_time);
    
    % find the begining of each trial using the raw signal
    [d,fs,dt,tc,t]=readhtk('ipad.htk');
    chan_ipad=d;
    plot(d)
    hold on
    THRESH =  input('thresh');
    START =  15% input('start');
    inds = find(chan_ipad(START:end-15)>=THRESH);
    [pos,n] = evFindGroups(inds,500,1);
    plot(inds(pos(1,1:end-1))+START,THRESH, '*r')
    
    ipad_ON_time = (inds(pos(1,1:end-1))+START)/(fs/1000);
    ipad_OFF_time = (inds(pos(1,2:end))+START)/(fs/1000);
    DELETE =  input('delete');
    GOOD= setdiff([1:length(ipad_ON_time)],DELETE)
    ecog(1).ipad_ON_time=int32(ipad_ON_time(GOOD));
    ecog(1).ipad_OFF_time= int32(ipad_OFF_time(GOOD));
    
else
%     ecog(1).rest_time=[];
%     ecog(1).active_time=[];
    ecog(1).ipad_ON_time=[];
    ecog(1).ipad_OFF_time=[];
    
end

%% save data
% use filename to create output filename
name=[filename '_ecog'];
filename=[filename,'_ecog.mat'];
save(filename,'name','ecog','M1_ch1','M1_ch2','aux','emg','Fs');

%% Add IPAD file info
name_ipad=[name(1:end-5) '_ipad.mat'];

try
    load(name_ipad);
    save(name ,'description','timestamp','-append');
end
