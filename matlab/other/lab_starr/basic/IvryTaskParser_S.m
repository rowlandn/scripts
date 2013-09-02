function EventsOut = IvryTaskParser_S(TimeStampVector, TaskChannel, SampleRate);

% This is the "bug-softening version" to be used with subject SCHW only.

% Function to turn the .abf data from the IvryTask into a list of
% events.  Events will be in data elements with respect to the initial
% data coming out of the .abf file; could be recoded into time if
% necessary.
%
% Written in March 2007 by John Schlerf, schlerf@berkeley.edu
  global TrialLengthThresh;
  global samplerate; % Making this global so that I can access it with helper functions
  
  %%%%%%%%%%%%%%%%%%%%%%
  % SampleRate is in samples/second, ie Hz (usually 20kHz)
  samplerate = SampleRate; % new variable name so I can pass it in
    
  %%%%%%
  % Constants that may eventually become parameters passed to the function
  TrialLengthThresh = 0.045;   % Minimum Trial Length, in Seconds
  VoltageThresh = 0.5;         % Threshold for a change in voltage
                               % necessary to trip our "change detector";
                               % noise is typically less than 0.05 upon
                               % initial inspection of the data, but
                               % raising this as high as 1 should still
                               % work just fine.
    
  % TaskChannel should move between 0 and +5 coming off the raw .abf file.  
  % However, ABFML spits it out as -1.2 to +51.  So, normalize it to keep
  % things in the right ballpark.  
  while max(TaskChannel) > 45
      TaskChannel = TaskChannel*0.1;
  end
  
  % The scheme is a little complex, so this routine should simplify it.
  
  
    
  % First thing, diff to find where it changes
  DiffedTask = diff(TaskChannel); 
  
  % Now let's use a quick thresholding procedure to ID
  % Rises and Falls:
  Rises = find(DiffedTask>VoltageThresh);
  Falls = find(DiffedTask<-VoltageThresh);
  % And be really anal in case the sampling rate catches 
  % this on a Voltage transition:
  Rises = Rises([1,find(diff(Rises)>1)+1]);
  Falls = Falls([1,find(diff(Falls)>1)+1]);
  
  % OK.  That lets me know where I should look for events, which is
  % good.  Now a trial should be a long stay in the "up" state, so I'll
  % do a quick subtraction to figure out when this happened.
  
  try % These should be the same size, so...
    UpTime = Falls-Rises;
  catch % But they might not be.  In that case, do the dumb for-loop way.
    warning(['The number of rises and falls is not equal; may be worth ' ...
          'double checking the output of this function']);
    for indx = 1:length(Rises);
      try
        UpTime(indx) = Falls(first(Falls > Rises(indx))) - Rises(indx);
      catch
        UpTime(indx) = -1;
      end
    end
  end
  
  
  %keyboard;
  
  % OK, now I know where I should see trials.  There should be more
  % "potential" trials than actual trials, as responses may flip this
  % into thinking theres three trials/trial.  But that's OK, we 
  % need to figure out what kind of trial they are anyway.
  
  % I can't think of a "smart" way to do this, so I'm going to use a few
  % crazy loops.  If there's a smarter way, I'd love to figure it out...

  RightDictionary = [1 10 100 11 12];
  LeftDictionary = [1 20 200 21 22];
  
  % Start by making empty arrays for all the variables I'll need:
  RightTrials = [];
  %RightStopTrials = [];
  LeftTrials = [];
  %LeftStopTrials = [];
  
  % These variables are going to be of length 5xN, as follows:
    % the first row holds the presentation of the fixation
    % the second row holds the presentation of the white stimulus
    % the third row holds the time it turned red 
      % (0 if it doesn't change color)
    % the fourth row holds the time the subject lifted off the button
      % (0 if they don't lift off)
    % the fifth row holds the time the subject pushed the button
      % (0 if they don't push)
  
  index = 1;
  PotentialTrials = lookForPotentialTrials(UpTime);
  KeepMe = PotentialTrials; % debugging thing
  while ~isempty(PotentialTrials)
      [LeftTrials, RightTrials, PotentialTrials] = parsetrial(Rises,Falls,UpTime,LeftTrials,RightTrials,PotentialTrials);
  end

  % OK, so I have all my timestamps, etc.  Right now time is in samples, so
  % I'll convert everything into the values in my TimeStamps vector that I
  % passed to this routine before making my output:
  home = GetTimes(TimeStampVector,LeftTrials(1,:)); home = [home GetTimes(TimeStampVector,RightTrials(1,:))];
  light = GetTimes(TimeStampVector,LeftTrials(2,:)); light = [light GetTimes(TimeStampVector,RightTrials(2,:))];
  stopon = GetTimes(TimeStampVector,LeftTrials(3,:)); stopon = [stopon GetTimes(TimeStampVector,RightTrials(3,:))];
  buttonup = GetTimes(TimeStampVector,LeftTrials(4,:)); buttonup = [buttonup GetTimes(TimeStampVector,RightTrials(4,:))];
  buttondown = GetTimes(TimeStampVector,LeftTrials(5,:)); buttondown = [buttondown GetTimes(TimeStampVector,RightTrials(5,:))];
  
  % The last thing I need is my trialtype parser, so let's do that here:
  type = ParseType(LeftTrials,{'A' 'C' 'E' 'Y'}); type = [type ParseType(RightTrials,{'B' 'D' 'F' 'Z'})];
  
  %%%%%%%%%%%%%%%%%%%%%%
  % JS - To name the fields differently, adjust below:
  
  EventsOut.home = home; % This is the start of fixation; at this time 
                         % subjects are known to be the "home" position.
                         % Ideally subjects should be in the home position
                         % during most of the down time, so this might need
                         % to be adjusted.  Also, zeros are inserted at 
                         % times, this does not indicate that they were not
                         % at home but that something was amiss.  Trials
                         % would not start if they were not at home.
  
  EventsOut.type = type; % This is my string of trial types.  Coded as
                         % follows:
                         % A - Left Go with Response
                         % B - Right Go with Response
                         % C - Left Stop, Successful
                         % D - Right Stop, Successful
                         % E - Left Stop, Failed
                         % F - Right Stop, Failed
                         % Y - Left Go with no response (Error)
                         % Z - Right Go with no response (Error)
                         
  EventsOut.cue = light; % This is when the cue came on
  
  EventsOut.stopsignal = stopon; % This is when the cue turned red.  By
                                 % necessity, this is 0 for trial types A 
                                 % and B.
                             
  EventsOut.lift = buttonup; % This is when the subject lifted off the 
                             % button.  By necessity, this is 0 for
                             % trial types C and D.
                                 
  EventsOut.press = buttondown; % This is when the subject pressed the 
                                % response button.  By necessity, this is
                                % 0 for trial types C and D, and may
                                % be 0 for types E and F.

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function typestring = ParseType(TrialInput,Dictionary)

% This will figure out what are go trials, and failed and successful stop
% trials

typestring = '';
for i = 1:size(TrialInput,2)
    if TrialInput(3,i) == 0
        if (TrialInput(4,i) ~= 0) & (TrialInput(5,i)~=0)
            typestring = [typestring Dictionary{1}];
        else
            typestring = [typestring Dictionary{4}];
        end
    elseif TrialInput(4,i) == 0
        typestring = [typestring Dictionary{2}];
    else
        typestring = [typestring Dictionary{3}];
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Events = GetTimes(TimeStampVector,Events)

% This translates event indices into TimeStamps; will leave zeros as they
% are.  Zeros indicate that no event was found for a given trial.

Events(find(Events)) = TimeStampVector(Events(find(Events)));

%%%%%%%%%%%%
% JS - This will change zeros to -1
% Events(find(Events == 0)) = -1; % or whatever you like

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PotentialTrials = lookForPotentialTrials(UpTimes)

% Simple behavior

global TrialLengthThresh;
global samplerate;

PotentialTrials = find(UpTimes > TrialLengthThresh * samplerate);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LeftTrials,RightTrials,PotentialTrials] = ...
    parsetrial(Rises,Falls,UpTime,LeftTrials,RightTrials,PotentialTrials)

    % the first row holds the presentation of the fixation
    % the second row holds the presentation of the white stimulus
    % the third row holds the time it turned red 
      % (0 if it doesn't change color)
    % the fourth row holds the time the subject lifted off the button
      % (0 if they don't lift off)
    % the fifth row holds the time the subject pushed the button
      % (0 if they don't push)
  
global TrialLengthThresh
global samplerate

% try
%     PotentialTrials = find(UpTime(index:index+6) > TrialLengthThresh*samplerate)+index-1;
% catch
%     PotentialTrials = find(UpTime(index:end) > TrialLengthThresh*samplerate)+index-1;
% end

FixTime = 0;
StartTime = 0;
StopSignalTime = 0;
ButtonUpTime = 0;
ButtonDownTime = 0;

% if isempty(PotentialTrials)
%     TrialRow = [0 0 0 0 0]';
% else
    StartIndx = PotentialTrials(1);
    
    Check = TrialType(UpTime(StartIndx-[3:-1:1]),Rises(StartIndx-[3:-1:1]));

    % This is "bug-softening" code.  It accomodates an apparent bug in
    % the Channel 10 code.
    dbFlag = 0;
    if isequal(Check,'MaybeRight')
        dbFlag = 1;
        
        try
            DownTimes = Rises(StartIndx+[1:2])-Falls(StartIndx+[0:1]);
            
            if (length(lookForPotentialTrials(UpTime(StartIndx+[0:2]))) == 3) ...
                & isempty(find(DownTimes > 0.005*samplerate))
                Check = 'Right';
            else
                Check = 'null';
            end
        catch
            Check = 'null';
        end
    end
    
    if (StartIndx > length(UpTime) - 1) % If we're at the end of the experiment... VERY special case...
        Rises(end+[1:3]) = Rises(end)+[5*UpTime(end) 10*UpTime(end) 20*UpTime(end)]; % make up 3 new ends
        UpTime(end+[1:3]) = [0.003 0.003 0.003]*samplerate; % make up 3 pulses
        Falls(end+[1:3]) = Rises(end-2:end)+UpTime(end-2:end); % ...and set 3 fall times.
    end
        
    
    if ~isequal(Check,'null') % If we found a real trial
        % Set up a variable to slim down the PotentialTrials variable.  The
        % logic is "what's the next element of PotentialTrials that we'll 
        % want to consider once we've taken this trial into account?"  
        % It'll be at least 2 in order to throw out the current value.
        ThrowFrom = 2; 
        
        % This will let me re-use all this code and have one "eval" line at
        % the end to assign the trial properly:
        if isequal(Check,'Right'), VarToUse = 'RightTrials'; else, VarToUse = 'LeftTrials'; end
      
        StartTime = Rises(StartIndx);
        try
            DownTimes = Rises(StartIndx+[1:3])-Falls(StartIndx+[0:2]);
        catch
            PotentialTrials = PotentialTrials(ThrowFrom:end);
            return
        end
        %keyboard
        if (DownTimes(1) < 0.005*samplerate)
            ButtonUpTime = Falls(StartIndx);
            if UpTime(StartIndx+1) > (TrialLengthThresh * samplerate) % PotentialTrials(2)
                ThrowFrom = ThrowFrom+1;
            end
            if (DownTimes(2) < 0.005*samplerate)
                ButtonDownTime = Falls(StartIndx+1);
                if UpTime(StartIndx+2) > (TrialLengthThresh * samplerate)
                    ThrowFrom = ThrowFrom+1;
                end
            end
        else 
            StopSignalTime = Falls(StartIndx); % That means the subject probably got a stop; code it!
            if ~isPulse(UpTime(StartIndx + 1)) % If the next UpTime is a pulse, then the subject
                                               % probably didn't press anything; IE not a failed
                                               % stop.  Otherwise...
                ButtonUpTime = Rises(StartIndx+1); % This will be close; these down pulses are a few
                                                   % samples long, so less than one ms of error
                if UpTime(StartIndx+1) > (TrialLengthThresh * samplerate)
                    ThrowFrom = ThrowFrom+1;
                end
                if (DownTimes(2) < 0.005*samplerate)
                    ButtonDownTime = Falls(StartIndx+1);
                    if UpTime(StartIndx+2) > (TrialLengthThresh * samplerate)
                        ThrowFrom = ThrowFrom+1;
                    end
                end
            end
        end
        try
            FixTime = Rises(find(isPulse(UpTime(StartIndx+[-10:0])))+StartIndx-10);
        catch
            FixTime = Rises(find(isPulse(UpTime(1:StartIndx))));
        end
        
        if isempty(FixTime), FixTime = [0 0]; elseif length(FixTime < 2), FixTime = [0 FixTime]; end
        
        % Set it to the pulse at the beginning of the counter
        if dbFlag, FixTime = FixTime(end); else, FixTime = FixTime(end-1); end
        
        TrialRow = [FixTime StartTime StopSignalTime ButtonUpTime ButtonDownTime]';
        
        PotentialTrials = PotentialTrials(ThrowFrom:end);
        eval([VarToUse '(:,end+1) = TrialRow;']);
    else
        PotentialTrials = PotentialTrials(2:end);
    end

return
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function type = TrialType(Pulses,Onsets)

% should return Right, Left, null, or MaybeRight

global samplerate

% type = 'valid'; % Placeholder variable

%OnsetDiff = Onsets(2)-Onsets(1);

if isempty(Pulses)
    keyboard
end

if ~(sum(isPulse(Pulses)) > 0) 
    if (Onsets(1) < 5000000) 
        type = 'null'; % Quick check
        %keyboard
        return;
    else
        % More bug softening.
        if isUpMarker(Pulses(3))
            if isUpMarker(Pulses(2))
                type = 'Left';
            else
                if Onsets(3)-(Onsets(2)+Pulses(2)) > 0.02*samplerate
                    type = 'Right';
                else
                    type = 'Left';
                end
            end
        else
            type = 'MaybeRight'; % There seems to be a bug in my Channel 10 code, in that there's not necessarily a delay after the counter.  So any code that depends on a 5 ms pulse being present will fail after trial # 64.
            return
        end
    end     
else
    if (isPulse(Pulses(3)))
        type = 'Right';
    elseif (isUpMarker(Pulses(2)) & isUpMarker(Pulses(3)))
        type = 'Left';
    elseif isPulse(Pulses(2)) & isUpMarker(Pulses(3))
        if Onsets(3) - Onsets(2) > 0.025*samplerate
            type = 'Right';
        else
            type = 'Left';
        end
    elseif (isPulse(Pulses(1)) & Onsets(1) > 5000000)
        type = 'MaybeRight';
    else
        % keyboard;
        type = 'null';
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boolout = isPulse(UpTimeLength)

global samplerate

boolout = (UpTimeLength < 0.010 * samplerate);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boolout = isUpMarker(UpTimeLength)

global samplerate

boolout = (UpTimeLength > 0.025 * samplerate)&(UpTimeLength < 0.035 * samplerate);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Index = first(BooleanVector)
    
    Index = find(BooleanVector);
    if ~isempty(Index)
        Index = Index(1);
    end
return
          