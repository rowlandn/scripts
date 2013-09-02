function events = ReadEventData( firstName, lastName, procedure, pass, track, event )
%
% Function that connects to the Guideline4000 database through ODBC and
% retrieves the data for the specified patient/procedure/pass/track/event
% If no event is provided, data for all events under a track is retreived
% If no track is provided, data for all events under a pass is retrieved
% 

% Parameters
driveOffset = 0;   % Normal use
%driveOffset = 29.950-9.130; % Drive has been reset at 9130. Initial distance to target was 29.950

% Switches
remote_db = 0;      % Set to 1 to connect to a remote database instead to the local one
display_rs_fields = 1;  % Switch for displaying the result set fields
normalize_waveforms = 0;  % When set to 1, a normalization of each waveform is performed; when 0, the actual relative amplitude is preserved.
scale = 4;  % The magnification factor for displying the waveforms
plot3d = 0;  % Set to 1 to do a 3D plot of the events and associated waveforms
plot_bkg = 1; % Switch for plotting the background RMS versus depth
wavelet_spike_removal = 0;
%wavelet_points = 65536; % In ADC samples
wavelet_points  = 4*65536; % In ADC samples
k_stddev_th     = 8.5;  % Number of standard deviations of the signal used as threshold for wavelet spike removal
k_snr_th        = 0.04;  % A correction factor for the wavelet spike removal threshold that weigths in the signal-to-noise ratio
max_wav_duration = 4; % Maximum waveform duration extracted from database, in secoonds. Set to 0 for unlimited waveform duration
plot_z = 1; % Switch for plotting the impedance versus depth

% Defaults
try
    firstName;
    lastName;
catch
    %error('At least first and last name must be provided.');
    
%     firstName = 'Dummy';
%     lastName  = 'Patient 2';
%     procedure = 'Procedure 1';
    %pass      = 'Pass 1 Right';    % 1 waveform
    %pass      = 'Pass 1 Left';      % 2 waveforms
    %track     = 'P';
    track     = '';
    %event     = 'Event 1';
    event     = '';

    firstName = 'Dummy';
    lastName  = 'Patient 3';
    procedure = 'Procedure 1';
    pass      = 'Pass 1';      % 2 waveforms
    track     = 'C';
    event     = 'Snapshot - 300.0 sec';

end;

try
    lastName;
catch
    lastName = '';
end;

try
    procedure;
catch
    procedure = '';
end;

try
    pass;
catch
    pass = '';
end;

try
    track;
catch
    track = '';
end;

try
    event;
catch
    event = '';
end;

% Constants
nchan   = 128;              % Default number of channels
maxwav  = 32768;
maxzwav = 16384;
maxswav = 16384;


% Initialization
%close all
figure(1);clf;
figure(2);clf;
figure(3);clf;
%figure(4);clf;
colordef(0,'black');

% Use java for retrieving information from database
%clear java
%javaclasspath

import java.sql.*
import java.io.*

%%%% Extract first patient in the list, its first procedure, first pass and first track from database
java.lang.Class.forName('sun.jdbc.odbc.JdbcOdbcDriver');
if (remote_db)
    url='jdbc:odbc:GuidelineDBRemote';
    con=DriverManager.getConnection(url, 'sa', 'saPWD');
    conWav=DriverManager.getConnection(url, 'sa', 'saPWD');
else
    url='jdbc:odbc:GuidelineDB';
    con=DriverManager.getConnection(url, '', '');
    conWav=DriverManager.getConnection(url, '', '');
end;
stmt=con.createStatement;
stmtWav=conWav.createStatement;
query=['SELECT Procedures.*, Passes.*, Tracks.TrackName, Tracks.Channel, Tracks.OffsetRL, Tracks.OffsetAP, Tracks.OffsetDV,' ...
       ' EvtID = CONVERT(nvarchar(36),Events.EventID),' ...
       ' Events.*, Events.TextNote AS EventNote,' ...
       ' ScanBasedMicroPointID = CONVERT(nvarchar(36),ScanBasedMicroLocations.PointID), ScanBasedMicroLocations.*,' ...
       ' WavID = CONVERT(nvarchar(36),Events.WaveformID),' ...
       ' AuxWavID = CONVERT(nvarchar(36),Events.AuxiliaryWaveformID),' ...
       ' ZWavID = CONVERT(nvarchar(36),ImpedanceMeasurements.WaveformID),' ...
       ' StimWavID = CONVERT(nvarchar(36),Stimulations.WaveformID),' ...
       ' Stimulations.*,' ...
       ' ImpedanceMeasurements.*' ...
       ' FROM Events' ...
       ' INNER JOIN Points AS ScanBasedMicroLocations ON Events.ScanBasedMicroLocationID = ScanBasedMicroLocations.PointID' ...
       ' LEFT OUTER JOIN ImpedanceMeasurements ON Events.ImpedanceMeasurementID = ImpedanceMeasurements.ImpedanceMeasurementID' ...
       ' LEFT OUTER JOIN Stimulations ON Events.StimulationID = Stimulations.StimulationID' ...
       ' INNER JOIN Tracks ON Events.TrackID=Tracks.TrackID' ...
       ' INNER JOIN Points ScanBasedTrackTops ON (Tracks.ScanBasedTrackTopLocationID=ScanBasedTrackTops.PointID)' ...
       ' INNER JOIN Passes ON Tracks.PassID=Passes.PassID INNER JOIN Procedures ON Passes.ProcedureID=Procedures.ProcedureID' ...
       ' INNER JOIN Patients ON Procedures.PatientID=Patients.PatientID' ...
       ' WHERE Patients.FirstName=''' firstName ''' AND Patients.LastName=''' lastName ''''];
       if (~isempty(procedure))
            query = [query ' AND Procedures.ProcedureName=''' procedure ''''];
       end;
       if (~isempty(pass))
            query = [query ' AND Passes.PassName=''' pass  ''''];
       end;
       if (~isempty(track))
            query = [query ' AND Tracks.TrackName=''' track ''''];
       end;
       if (~isempty(event))
            query = [query ' AND Events.Caption=''' event ''''];
       end;
       %query = [query ' ORDER BY Procedures.ProcedureName, Passes.PassName, Tracks.Channel, Events.DrivePosition DESC'];    % Use this for distance traveled
       query = [query ' ORDER BY Procedures.ProcedureName, Passes.PassName, Tracks.Channel, Events.DrivePosition'];    % Use this for distance to target


disp(query);
rs=stmt.executeQuery(query);
events=[];
eventNum = 0;
wavNum = 0;
zWavNum = 0;
stimWavNum = 0;
while (rs.next)
    eventNum = eventNum + 1;
    % List column names
    if (display_rs_fields)
        md = rs.getMetaData;
        for i=1:md.getColumnCount
            disp(sprintf('  %.2d: %s',i,char(md.getColumnName(i))));
            %disp(md.getColumnType(i));
        end;
        display_rs_fields = 0;
    end;

    % Must read data in the order in which it shows in the result set (see
    % output from above for loop)
    %patientLastName = char(rs.getString('LastName'))
    procedureName = char(rs.getString('ProcedureName'));
    passName = char(rs.getString('PassName'));
    trackName = char(rs.getString('TrackName'));
    channel = double(rs.getInt('Channel'));
    offsetRL    = double(rs.getDouble('OffsetRL'));
    offsetAP    = double(rs.getDouble('OffsetAP'));
    offsetDV    = double(rs.getDouble('OffsetDV'));
    eventID     = char(rs.getString('EvtID'));
    eventName   = char(rs.getString('Caption'));
    microDistanceToTarget = -double(rs.getDouble('MicroDistanceToTarget')) + driveOffset;
    macroDistanceToTarget = -double(rs.getDouble('MacroDistanceToTarget')) + driveOffset;
    isMacroEvent = double(rs.getInt('IsMacroEvent'));
    if (isMacroEvent==1)
        distanceToTarget = macroDistanceToTarget;
    else
        distanceToTarget = microDistanceToTarget;
    end;
    eventNote = char(rs.getString('EventNote'));
    %size(eventNote)
    %strtok(eventNote)
    pointID  = char(rs.getString('ScanBasedMicroPointID'));
    RL    = double(rs.getDouble('RL'));
    AP    = double(rs.getDouble('AP'));
    DV    = double(rs.getDouble('DV'));
    waveformID          = char(rs.getString('WavID'));
    auxWaveformID       = char(rs.getString('AuxWavID'));
    zWaveformID         = char(rs.getString('ZWavID'));
    stimWaveformID      = char(rs.getString('StimWavID'));
    stimulationCurrent  = double(rs.getDouble('StimulationCurrent'));
    Z                   = double(rs.getDouble('Impedance'));
    disp(sprintf('%d :\t%s\t%s\tch:%d(%s)\tI=%.3f\ttarget:%.3f\tRL:%.3f\tAP:%.3f\tDV:%.3f (spike: %s/ Z: %s / Stim: %s)',eventNum,eventName,passName,channel,trackName,stimulationCurrent,distanceToTarget,RL,AP,DV,waveformID,zWaveformID,stimWaveformID));
    
    events(eventNum).procedureName = procedureName;
    events(eventNum).passName = passName;
    events(eventNum).trackName = trackName;
    events(eventNum).channel = channel;
    events(eventNum).eventName = eventName;
    events(eventNum).microDistanceToTarget = microDistanceToTarget;
    events(eventNum).macroDistanceToTarget = macroDistanceToTarget;
    events(eventNum).RL = RL;
    events(eventNum).AP = AP;
    events(eventNum).DV = DV;
    events(eventNum).Note = eventNote;
    events(eventNum).Z = Z;
    
    % Recorded spike waveform
    % Retreiving the waveform from the database takes quite a bit. Rather
    % than retreiving the waveform everytime, we attempt to load previously
    % saved (or 'cached') waveforms from disk;
    wavD2T = [];
    if (~isempty(waveformID) && ~strcmp(waveformID,'00000000-0000-0000-0000-000000000000'))
        wavNum = wavNum + 1;
        wavD2T(wavNum)=distanceToTarget;
        %wavfilename = sprintf('%s.wav',waveformID);    % Use the waveform ID
        wavfilename = sprintf('%s %s %s %s %s %s D%.2f.wav',firstName,lastName,procedureName,passName,trackName,eventName,microDistanceToTarget);
        matfilename = sprintf('%s %s %s %s %s %s D%.2f.mat',firstName,lastName,procedureName,passName,trackName,eventName,microDistanceToTarget);
        auxmatfilename = sprintf('%s %s %s %s %s %s D%.2f.aux.mat',firstName,lastName,procedureName,passName,trackName,eventName,microDistanceToTarget);
        if (exist(matfilename,'file')==2)
            % Waveform exists, load it from disk
            load(matfilename);
        else
            % Retrieve waveform (an APM blob) from database
            queryWav=['SELECT Waveforms.WaveformData, Waveforms.WaveformURI FROM Waveforms WHERE Waveforms.WaveformID=''' waveformID ''''];
            disp(queryWav);
            rsWav=stmtWav.executeQuery(queryWav);
            if (rsWav.next)
                wb=rsWav.getBytes('WaveformData');
                waveformURI = char(rsWav.getString('WaveformURI'));
                if (~isempty(wb))
                    apmfilename = sprintf('%s %s %s %s %s %s D%.2f.apm',firstName,lastName,procedureName,passName,trackName,eventName,microDistanceToTarget);
                    fapm=fopen(apmfilename,'wb');
                    fwrite(fapm,wb,'int8');
                    fclose(fapm);
                else
                    apmfilename = waveformURI;
                end;
                t=APMReadData(apmfilename);
                % Save waveform as a mat file
                if (~isempty(t))
                    save(matfilename,'t');
                    waveform = t.channels(channel).continuous;
                    sampfreq(channel) = t.channels(channel).sampling_frequency;
                    if (~isempty(waveform))
                        if (normalize_waveforms)
                            % The voltage calibration must be adjusted to compensate for
                            % the normalization
                            % ...
                            waveform = waveform / max(abs(waveform));
                        else
                            waveform = waveform / maxwav;
                        end;
                        wavwrite(waveform,sampfreq(channel),16,wavfilename);
                    end;
                end;
            end;
        end;
        if (exist(auxmatfilename,'file')==2)
            % Waveform exists, load it from disk
            load(auxmatfilename);
        else
            % Retrieve auxiliary waveform (an APM blob) from database
            queryWav=['SELECT Waveforms.WaveformData, Waveforms.WaveformURI FROM Waveforms WHERE Waveforms.WaveformID=''' auxWaveformID ''''];
            disp(queryWav);
            rsWav=stmtWav.executeQuery(queryWav);
            if (rsWav.next)
                wb=rsWav.getBytes('WaveformData');
                %size(wb);
                waveformURI = char(rsWav.getString('WaveformURI'));
                if (~isempty(wb))
                    auxapmfilename = sprintf('%s %s %s %s %s %s D%.2f.aux.apm',firstName,lastName,procedureName,passName,trackName,eventName,microDistanceToTarget);
                    fapm=fopen(auxapmfilename,'wb');
                    fwrite(fapm,wb,'int8');
                    fclose(fapm);
                else
                    auxapmfilename = waveformURI;
                end;
                taux=APMReadData(auxapmfilename);
                save(auxmatfilename,'taux');
            end;
        end;
        % Detect spikes, based on the parameters stored in the apm file
        % First check if discrim params seem to have been set correctly
        maxw=max(abs(t.channels(channel).continuous));
        if ( (abs(t.channels(channel).trig_level)<maxw) && (abs(t.channels(channel).win_low)<maxw) )
        else
            % Choose some defaults
            if (abs(max(t.channels(channel).continuous))>abs(min(t.channels(channel).continuous)))
                t.channels(channel).trig_level = 0.3 * maxw;
                t.channels(channel).win_low = 0.6 * maxw;
                t.channels(channel).win_high = 1.1 * maxw;
                t.channels(channel).trig_slope = 1;
            else
                t.channels(channel).trig_level = - 0.3 * maxw;
                t.channels(channel).win_low = - 0.6 * maxw;
                t.channels(channel).win_high = - 1.1 * maxw;
                t.channels(channel).trig_slope = -1;
            end;
        end;
        if (abs(t.channels(channel).win_high)>1.12*maxw)
            t.channels(channel).win_high=1.12*t.channels(channel).win_high/maxw;
        end;
        %eventNote;
        %t.channels(channel).win_low
        %t.channels(channel).win_high
        if (~isempty(strfind(lower(eventNote),'sorted')))
            % Sort spikes
            t=WDGL4K(t,channel); % Use level discriminator
        else
            t.channels(channel).timestamps=[];
        end;
        
        mfr=length(t.channels(channel).timestamps)/(length(t.channels(channel).continuous)/t.channels(channel).sampling_frequency);
        disp(sprintf('Mean firing rate over the selected interval is %.2f sp/sec',mfr));
        
        waveform = t.channels(channel).continuous/maxwav;
        sampfreq(channel) = t.channels(channel).sampling_frequency;
        events(eventNum).waveform = waveform;   % The waveform is scaled from 0 to 1 (as this is required by wavwrite)
        events(eventNum).samplingFrequency = sampfreq(channel); 
        
        % Calculate the mean power of the signal (matching the squared standard deviation), including spikes.
        % A further step includes removing the spikes and calculating the background activity
        % Perform any additional filtering here
        if (wavelet_spike_removal)
            Jmin = 2; 
            options.wavelet_type = 'daubechies';
            options.wavelet_vm = 4;
            options.ti = 0; % Set this to 1 to test for translation-invariant wavelet
            % wavelet transform
            nwp=wavelet_points;
            while nwp>length(waveform)
                nwp=nwp/2;
            end;
            %y = detrend(waveform,'constant');   % Substract DC offset, if any (due to improprer calibration or variable gain offset)
            y=waveform(1:nwp);
            y = detrend(y,'constant');   % Substract DC offset, if any (due to improprer calibration or variable gain offset)
            %wthresh=nstddevT*std(y,1);  % nstddevT standard deviations 
            wthresh=k_stddev_th * (1 - k_snr_th * max(abs(y)) / std(y,1)) * std(y,1);  % nstddevT standard deviations corrected for high SNR ratio
            %[max(abs(waveform)) std(waveform,1)]
            yw = perform_wavelet_transform(y,Jmin,+1, options); 
            % perform thresholding
            ywT = yw .* (abs(yw)>wthresh);
            % inverse transform
            yT = perform_wavelet_transform(ywT,Jmin,-1, options); 
            ynoise = y - yT; 
            stddev2=std(ynoise,1)^2;
        else
            stddev2=std(waveform,1)^2;
        end;
        events(eventNum).stddev2 = stddev2; 
        
        
        figure(1);
        set(1,'Color','k');
        subplot('Position',[0.2 0.1 0.8 0.8]);
        hold on;
        title(sprintf('Recorded Waveforms Track %s',trackName),'Visible','On');
        set(gca,'Visible','Off');
        % Plot waveform
        plot(waveform*scale-2*wavNum,'g','Visible','On');  % Display each waveform, shifted by its index
        if (wavelet_spike_removal)
            plot(ynoise*scale-2*wavNum,'-','Color',[1 0.4 0.1],'Visible','On');  % Display each waveform, shifted by its index
        end;
        % Plot timestamps
        for i=1:length(t.channels(channel).timestamps)
            ts=t.channels(channel).timestamps(i)-t.channels(channel).start_continuous+1;
            plot([ts ts],[0.8 1]-2*wavNum,'r-','LineWidth',1);
        end;
        plot([1 length(waveform)],[0 0],'-','Color',[0.7 0.7 0.7]); % Zero line
        plot([1 length(waveform)],[t.channels(channel).trig_level/maxwav t.channels(channel).trig_level/maxwav]*scale-2*wavNum,'-','Color',[0.7 0.4 0.4]);
        plot([1 length(waveform)],[t.channels(channel).win_low/maxwav t.channels(channel).win_low/maxwav]*scale-2*wavNum,'-','Color',[0.7 0.7 0.7]);
        plot([1 length(waveform)],[t.channels(channel).win_high/maxwav t.channels(channel).win_high/maxwav]*scale-2*wavNum,'-','Color',[0.7 0.7 0.7]);
        text(1,-2*wavNum+0.5,[waveformID '.wav'],'FontSize',5,'Color',[0.5 0.5 0.5]);
        text(1,-2*wavNum+0.7,sprintf('%s   ch:%d    %.3f  MFR=%.0f sp/sec  {\\sigma}^2=%.1f',eventName,channel,microDistanceToTarget,mfr,1e6*stddev2),'Color','w','FontSize',8);
        hold off;
        subplot('Position',[0.05 0.1 0.1 0.8]);
        hold on;
        box on;
        set(gca,'Visible','Off');
        set(gca,'FontSize',7);
        plot(channel,-distanceToTarget,'r.','MarkerSize',20);
        text(0.42+channel,-distanceToTarget,sprintf('%.3f',distanceToTarget),'FontSize',7);
        %text(channel,-distanceToTarget,sprintf('%.3f',trackName),'FontSize',7);
        hold off
        drawnow;
        pause(0.1);
    end;
    
    % Recorded impedance measurements waveform
    % Retreiving the waveform from the database takes quite a bit. Rather
    % than retreiving the waveform everytime, we attempt to load previously
    % saved (or 'cached') waveforms from disk;
    if (~isempty(zWaveformID))
        zWavNum = zWavNum + 1;
        %zWavFilename = sprintf('%s.wav',zWaveformID);    % Use the waveform ID
        zWavFilename = sprintf('%s %s %s %s %s %s D%.2f Z.wav',firstName,lastName,procedureName,passName,trackName,eventName,distanceToTarget);
        if (exist(zWavFilename,'file')==2)
            % Waveform exists, load it from disk
            zWaveform=wavread(zWavFilename);
        else
            zWaveform=[];
            % Retrieve waveform (an APM blob) from database
            queryWav=['SELECT Waveforms.WaveformData FROM Waveforms WHERE Waveforms.WaveformID=''' zWaveformID ''''];
            disp(queryWav);
            rsWav=stmtWav.executeQuery(queryWav);
            if (rsWav.next)

                wd=rsWav.getBinaryStream('WaveformData');

                wds=java.io.DataInputStream(wd);   % This works for big endian representation
                %wds=com.mindprod.ledatastream.LEDataInputStream(wd); % This works for little endian (Intel) representation
                % Data comes as a vector of integers representing the SC ADC
                % values
                while (1)
                    try
                        k=wds.readInt;
                        zWaveform = [zWaveform k];
                    catch
                        disp('Done processing impedance measurement waveform');
                        break;
                    end;
                end;
                wds.close;
                wd.close;
                if (~isempty(zWaveform))
                    if (normalize_waveforms)
                        zWaveform = zWaveform / max(abs(zWaveform));
                        % The voltage calibration must be adjusted to compensate for
                        % the normalization
                        % ...
                    else
                        zWaveform = zWaveform / maxzwav;
                    end;
                    wavwrite(zWaveform,40000,16,zWavFilename);
                end;
            else
                disp('Impedance measurement waveform not found in database.');
            end;
        end;
        if (~isempty(zWaveform))
            figure(2);
            set(2,'Color','k');
            hold on;
            title('Impedance Measurement Waveform');
            set(gca,'Visible','Off');
            plot(zWaveform+2*zWavNum,'g','Visible','On');  % Display each waveform, shifted by its index
            text(1,2*zWavNum+0.5,[zWaveformID '.wav'],'FontSize',5,'Color',[0.5 0.5 0.5]);
            text(1,2*zWavNum+0.7,sprintf('%s   ch:%d    %.3f',eventName,channel,microDistanceToTarget),'Color','w','FontSize',8);
            hold off;
            drawnow;
            pause(0.1);
        end;
    end;
    
    % Recorded impedance measurements waveform
    % Retreiving the waveform from the database takes quite a bit. Rather
    % than retreiving the waveform everytime, we attempt to load previously
    % saved (or 'cached') waveforms from disk;
    if (~isempty(stimWaveformID))
        stimWavNum = stimWavNum + 1;
        %stimWavFilename = sprintf('%s.wav',stimWaveformID);    % Use the waveform ID
        stimWavFilename = sprintf('%s %s %s %s %s %s D%.2f Stim.wav',firstName,lastName,procedureName,passName,trackName,eventName,distanceToTarget);
        if (exist(stimWavFilename,'file')==2)
            % Waveform exists, load it from disk
            stimWaveform=wavread(stimWavFilename);
        else
            stimWaveform=[];
            % Retrieve waveform (an APM blob) from database
            queryWav=['SELECT Waveforms.WaveformData FROM Waveforms WHERE Waveforms.WaveformID=''' stimWaveformID ''''];
            disp(queryWav);
            rsWav=stmtWav.executeQuery(queryWav);
            if (rsWav.next)

                wd=rsWav.getBinaryStream('WaveformData');

                wds=java.io.DataInputStream(wd);   % This works for big endian representation
                %wds=com.mindprod.ledatastream.LEDataInputStream(wd); % This works for little endian (Intel) representation
                % Data comes as a vector of integers representing the SC ADC
                % values
                while (1)
                    try
                        k=wds.readInt;
                        stimWaveform = [stimWaveform k];
                    catch
                        disp('Done processing stimulation waveform');
                        break;
                    end;
                end;
                wds.close;
                wd.close;
                if (~isempty(stimWaveform))
                    if (normalize_waveforms)
                        % The voltage calibration must be adjusted to compensate for
                        % the normalization
                        % ...
                        stimWaveform = stimWaveform / max(abs(stimWaveform));
                    else
                        stimWaveform = stimWaveform / maxswav;
                    end;
                    wavwrite(stimWaveform,40000,16,stimWavFilename);
                end;
            else
                disp('Stimulation waveform not found in database.');
            end;
        end;
        if (~isempty(stimWaveform))
            figure(3);
            set(3,'Color','k');
            hold on;
            title('Stimulation Waveform');
            set(gca,'Visible','Off');
            plot(stimWaveform+2*stimWavNum,'g','Visible','On');  % Display each waveform, shifted by its index
            text(1,2*stimWavNum+0.5,[stimWaveformID '.wav'],'FontSize',5,'Color',[0.5 0.5 0.5]);
            text(1,2*stimWavNum+0.7,sprintf('%s   ch:%d    %.3f',eventName,channel,microDistanceToTarget),'Color','w','FontSize',8);
            hold off;
            drawnow;
            pause(0.1);
        end;
    end;
    
    % Now plot stuff in 3D
    if (plot3d)
        figure(4);
        set(4,'Color','k');
        hold on;
        %set(gca,'CameraPosition',[-61.5585 323.736 -118.479]); 
        %set(gca,'CameraPosition',[226.555 -193.673 -120.077]);
        set(gca,'CameraPosition',[15 346.41 0]);
        xlabel('RL');
        ylabel('AP');
        zlabel('DV');
        plot3(RL,AP,DV,'r.','MarkerSize',10);
        % Plot the MRI slice, if present
        if (~isempty(who('D')))
            scansize = 300;     % The scan size, in mm
            %surf(linspace(-0.5*scansize,0.5*scansize,size(D,1)),linspace(-0.5*scansize,0.5*scansize,size(D,2))
        end;
        if (~isempty(waveformID))
            % Show the waveform thumbnail just next to the event (2 mm laterally)
            mmpersec=2; % Number of mm per second in the waveform thumbnails
            mmperwav=16; % Number of mm of the thumbnail height
            h=plot3(RL+2+mmpersec*linspace(0,length(waveform)/sampfreq(channel),length(waveform)),AP+zeros(1,length(waveform)),DV+mmperwav*waveform,'g-');
        end;
        hold off;
        drawnow;
        pause(0.1);
    end;
    
end;

con.close;
disp(sprintf('Processed %d events',eventNum));

% Plot the squared standard deviation (matching the RMS value of the signal)
if (plot_bkg)
    d2t=[];
    stddev2=[];
    for i=1:length(events)
        if ~isempty(events(i).stddev2)
            d2t=[d2t;events(i).microDistanceToTarget];
            stddev2=[stddev2;events(i).stddev2];
        end;
    end;
    stddev2=stddev2*1e6;
    figure(4);clf;
    hold on;
    title(sprintf('Background Signal Power %s %s %s %s Track %s',firstName,lastName,procedure,pass,track));
    xlabel('Distance to target (mm)');
    ylabel('Background Power (a.u.)');
    plot([0 0],[0 max(stddev2)],'w--');
    plot(d2t,stddev2,'go-');
    hold off;
end;

% Plot the impedance vs depth
if (plot_z)
    d2t=[];
    Z=[];
    for i=1:length(events)
        if (events(i).Z > 0)
            d2t=[d2t;events(i).microDistanceToTarget];
            Z=[Z;events(i).Z];
        end;
    end;
    figure(5);clf;
    hold on;
    title(sprintf('Impedance %s %s %s %s Track %s',firstName,lastName,procedure,pass,track));
    xlabel('Distance to target (mm)');
    ylabel('Impedance (M{\Omega})','Interpreter','TeX');
    plot([0 0],[0 max(Z)],'w--');
    plot(d2t,Z,'go-');
    hold off;
end;

