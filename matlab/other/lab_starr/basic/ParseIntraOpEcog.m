ParseIntraOpEcog( firstName, lastName, procedure, pass, track, event )

% read event data from Guideline4000 database
try
    events = ReadEventData( firstName, lastName, procedure, pass, track, event );
    t1 = events.waveform;
catch
    % empty aux channels give errors in ReadEventData
    % use continuous channel data contained in mat until error is fixed
    filename = [firstName ' ' lastName ' '...
        procedure ' ' pass ' '...
        track ' ' event ' '];
    
end
