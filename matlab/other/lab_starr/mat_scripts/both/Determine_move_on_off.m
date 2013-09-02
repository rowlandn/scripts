function Determine_move_on_off(name)

% Determine movement onset and offset of movement using emg/accel/task button timing data.

% Created by Cora (3/15/2012)

%% load data

load(name)

%create time matching time vector
    nsamples = length(aux.chan(2).raw); %#ok<NODEF>
    T = 1/(Fs);
    time = 0:T:T*(nsamples-1);
   %[move_onset,move_offset,bad_move_onset,bad_move_offset] =  DetectMove_ON_OFF(time,aux.chan(2).raw,aux.chan(3).raw,emg.chan(1).raw,ecog.go_time/Fs,ecog.ipad_OFF_time/Fs);
 [move_onset,move_offset,bad_move_onset,bad_move_offset] =  DetectMove_ON_OFF(time,aux.chan(2).raw,ecog.go_time/Fs,ecog.ipad_OFF_time/Fs);
 
ecog.active_time = int32(move_onset*Fs);
ecog.rest_time = int32(move_offset*Fs);
ecog.bad_move_on_time = int32(bad_move_onset*Fs);
ecog.bad_move_off_time = int32(bad_move_offset*Fs);
save(name,'ecog','-append')

