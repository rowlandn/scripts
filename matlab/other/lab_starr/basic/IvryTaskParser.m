function EventsOut = IvryTaskParser(TaskChannel);

%

% Function to turn the .abf data from the IvryTask into a list of
% events.  Events will be in data elements with respect to the initial
% data coming out of the .abf file; could be recoded into time if
% necessary.
%
% Written in March 2007 by John Schlerf, schlerf@berkeley.edu
  global TrialLengthThresh;
  
  %%%%%%
  % Constants that may eventually become parameters passed to the function
  SampleRate = 20000;           % samples/second
  TrialLengthThresh = 0.150;   % Minimum Trial Length, in Seconds
  VoltageThresh = 0.1;         % Threshold for a change in voltage
                               % necessary to trip our "change detector";
                               % noise is typically less than 0.05 upon
                               % initial inspection of the data, but
                               % raising this as high as 1 should still
                               % work just fine.
    
  global samplerate;                             
  samplerate = SampleRate;  % Stupid hack, but functional.
  % TaskChannel should move between 0 and +5 coming off
  % the raw .abf file.  The scheme is a little complex, so 
  % this routine should simplify it.
    
  % First thing, diff to find where it changes
  DiffedTask = diff(TaskChannel); 
  
  % Now let's use a quick thresholding procedure to ID
  % Rises and Falls:
  Rises = find(DiffedTask>VoltageThresh);
  Falls = find(DiffedTask<-VoltageThresh);
  % And be really anal in case the sampling rate catches 
  % this on a Voltage transition:
  Rises = Rises([1;find(diff(Rises)>1)+1]);
  Falls = Falls([1;find(diff(Falls)>1)+1]);
  
  % OK.  That lets me know where I should look for events, which is
  % good.  Now a trial should be a long stay in the "up" state, so I'll
  % do a quick subtraction to figure out when this happened.
  
  try % These should be the same size, so...
    UpTime = Falls-Rises;
  catch % But they might not be.  In that case, do the dumb for-loop way.
    disp(['The number of rises and falls is not equal; may be worth ' ...
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
  
  for index = 1:length(UpTime)-2
      switch(TrialType(UpTime(index:index+2),Rises(index:index+2)))
          
          case 'Right'
              RightTrials(:,end+1) = parsetrial(Rises,Falls,UpTime,index);
          %case 'RightStop'
              %RightStopTrials(:,end+1) = parsetrial(Rises,Falls,UpTime,index);
          case 'Left'
              LeftTrials(:,end+1) = parsetrial(Rises,Falls,UpTime,index);
          %case 'LeftStop'
              %LeftStopTrials(:,end+1) = parsetrial(Rises,Falls,UpTime,index);
      end
  end

  keyboard;

  EventsOut = zeros(size(TaskChannel));
  
  EventsOut = AddOnsets(EventsOut,LeftTrials,LeftDictionary);
  EventsOut = AddOnsets(EventsOut,RightTrials,RightDictionary);
  
  keyboard;
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatrixOut = AddOnsets(MatrixOut,Trials,Dictionary);

    for Row = 1:size(Trials,1);
        Indices = Trials(Row,:);
        Indices = Indices(find(Indices));
        MatrixOut(Indices) = Dictionary(Row);
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TrialRow = parsetrial(Rises,Falls,UpTime,index)
  % Returns a vector of length length 5, as follows:
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

try
    PotentialTrials = find(UpTime(index:index+6) > TrialLengthThresh*samplerate)+index-1;
catch
    PotentialTrials = find(UpTime(index:end) > TrialLengthThresh*samplerate)+index-1;
end

FixTime = 0;
StartTime = 0;
StopSignalTime = 0;
ButtonUpTime = 0;
ButtonDownTime = 0;

if isempty(PotentialTrials)
    TrialRow = [0 0 0 0 0]';
else
    StartIndx = PotentialTrials(first(PotentialTrials));
    StartTime = Rises(StartIndx);
    DownTimes = Rises(StartIndx+[1:3])-Falls(StartIndx+[0:2]);
    %keyboard
    if (DownTimes(1) < 0.005*samplerate)
        ButtonUpTime = Falls(StartIndx);
        if (DownTimes(2) < 0.005*samplerate)
            ButtonDownTime = Falls(StartIndx+1);
        end
    else %if ~isPulse(UpTime(StartIndx+1)) % So if the channel went down before they pressed a button and it's not followed by a pulse...
        
        StopSignalTime = Falls(StartIndx); % That means the subject probably got a stop; code it!
        if UpTime(StartIndx + 1) > 0.010*samplerate
            ButtonUpTime = Rises(StartIndx+1)-floor(0.002*samplerate);
            if (DownTimes(2) < 0.005*samplerate)
                ButtonDownTime = Falls(StartIndx+1);
            end
        end
    end
    TrialRow = [FixTime StartTime StopSignalTime ButtonUpTime ButtonDownTime]';
end

return
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function type = TrialType(Pulses,Onsets)

% should return RightGo, RightStop, LeftGo, LeftStop, or null

global samplerate

type = 'valid'; % Placeholder variable

OnsetDiff = Onsets(2)-Onsets(1);

if ~isPulse(Pulses(1))
    type = 'null';
else
    if (~isUpMarker(Pulses(2))) & (OnsetDiff > 0.15*samplerate)
        type = 'Right';
    elseif (isUpMarker(Pulses(2)) & isUpMarker(Pulses(3)))
        type = 'Left';
    elseif isUpMarker(Pulses(2)) & (Onsets(3)-Onsets(2) > 0.1*samplerate)
        if OnsetDiff > 0.025*samplerate
            type = 'Right';
        else
            type = 'Left';
        end
    else
        if ~isPulse(Pulses(3))
            %keyboard;
        end
        type = 'null';
    end
end

%     elseif (isUpMarker(Pulses(2)) & ((Onsets(2)-O
%     junk = diff(Onsets);
%     if (junk(1) > .007 * samplerate) & (junk(1) < .015 * samplerate)
%         if (Pulses(2) > 0.025*samplerate)&(Pulses(2)<0.035*samplerate)
%             First = 'Right'
%             EvalMe = sum(junk);
%         else
%             type = 'null';
%             EvalMe = 0;
%         end
%     else
%         First = 'Left';
%         EvalMe = junk(1);
%     end;
%     if (EvalMe > .045 * samplerate) & (EvalMe < .055 * samplerate)
%         PulseToLookAt = Pulses(3-isequal(First,'Left'));
%         if (PulseToLookAt > 0.025*samplerate)&(PulseToLookAt < 0.035*samplerate)
%             Second = 'Stop';
%         else
%             type = 'null';
%         end
%     else
%         Second = 'Go';
%     end;
%     if ~isequal(type,'null')
%         type = [First Second];
%     end
% end
% 
% if ~isequal(type,'null')
%     % one last double check
%     if isequal(First,'Right');
%         if (Pulses(2) < .055*samplerate) | (Pulses(2) > .065*samplerate)
%             type = 'null';
%         end
%         if isequal(Second,'Stop')
%             if (Pulses(3) < .055*samplerate) | (Pulses(3) > .065*samplerate)
%                 type = 'null';
%             end
%         end
%     else
%         if isequal(Second,'Stop')
%             if (Pulses(2) < 0.055*samplerate) | (Pulses(2) > .065*samplerate)
%                 type = 'null';
%             end
%         else
%             if Pulses(2) < .200*samplerate
%                 type = 'null';
%             end
%         end
%     end
% end

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
          