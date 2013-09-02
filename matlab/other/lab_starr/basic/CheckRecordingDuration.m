
maxDuration = 3599;

% Use java for retrieving/updating information in the database
import java.sql.*
import java.io.*

java.lang.Class.forName('sun.jdbc.odbc.JdbcOdbcDriver');
url='jdbc:odbc:GuidelineDB';
con=DriverManager.getConnection(url, '', '');
%url='jdbc:odbc:GuidelineDBRemote';                      % Use this for remote connections
%con=DriverManager.getConnection(url, 'sa', 'saPWD');    % Use this for remote connections
stmt=con.createStatement;

% Check the WaveformDuration field in the Events table
query='SELECT Events.* FROM Events WHERE WaveformDuration>3599';
%disp(query);
rs=stmt.executeQuery(query);
eventNum = 0;
waveformDuration = [];
while (rs.next)
    eventNum = eventNum + 1;
    waveformDuration(eventNum)  = double(rs.getDouble('WaveformDuration'));
end;
rs.close;
disp(sprintf('Found %d events with large waveform duration:',eventNum));
disp(waveformDuration.');
update = 'UPDATE Events SET WaveformDuration = 3599 WHERE WaveformDuration > 3599';
rs=stmt.executeUpdate(update);

% Check the Duration field in the Waveforms table
query='SELECT Waveforms.Duration FROM Waveforms WHERE Duration>3599';
%disp(query);
rs=stmt.executeQuery(query);
wavNum = 0;
waveformDuration = [];
while (rs.next)
    wavNum = wavNum + 1;
    waveformDuration(wavNum)  = double(rs.getDouble('Duration'));
end;
rs.close;
disp(sprintf('Found %d waveforms with large duration:',wavNum));
disp(waveformDuration.');
update = 'UPDATE Waveforms SET Duration = 3599 WHERE Duration > 3599';
rs=stmt.executeUpdate(update);


disp('Done.');
con.close;

disp(sprintf('Processed %d events and %d waveforms',eventNum,wavNum));


