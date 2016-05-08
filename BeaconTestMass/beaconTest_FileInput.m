% Matlab code for checking beacons.  Looking for 'Bad Beacons'
% This program accepts the serial input from a guidance light board that is
% outputting custom serial data that has the HID,RSSI,...
%
% Sections commented as required

% ***********Environment setup************
% close any currently open figures
close all
instrreset; % Close any open serial ports.
fclose('all');  % Close any open files
clear all
% open a figure to allow keyboard capture to end the infinite loop and to
% plot data in live view.
figHandle = figure(1);
%******************************************

% ******Define local variables***********
beaconMax = 30; % The maximum number of beacons expected in a packet
beconHID = zeros(1,beaconMax); % An array for the HID
beaconRSSI = zeros(1,beaconMax); % An array for the RSSI (note, will grow in a loop) 

maxRows = 1000; % The maximum number of rows that we expect to receive (note that if the actual number is larger then the system will auto grow the array

serialPortNumber = '4';   % The number of the target comm port
serialBaud = 115200;    % Target baud rate for the serial
serialInputBufferSize = 4048;   % defined size for the input buffer
serialFlowControl = 'hardware';  % defined flow control type

% Log file output settings.  Note that this will need to be changed if the
% program is ported to another computer.  If necessary, look to prompt the
% user for a file name
%logFilePath = 'C:\Users\Mathew\Documents\MATLAB\BeaconTestLogs\';
logFilePath = 'C:\Documents and Settings\Administrator\My Documents\MATLAB\Beacon testing\BeaconTestLogOutput\';
logFileName = strcat(timeDateNow,' beaconTestLog.csv');
logFileFullPath = strcat(logFilePath,logFileName);

escape = false; % Escape from infinite while loop tracker

runningAvWindow = 5;    % adjustable value for the running average.
%****Init variables to store measurements and details of the beacons. 
% Note that this has preallocated size but we will only work on the filled 
% entries.

beaconHID = zeros(beaconMax);   %Place to store the HIDs of the beacons.  zeroCorrectionCount = zeros(beaconMax); %tracks the number of times each beacon has it's value averaged

legendHID = cell(beaconMax);    % To store the names of HIDs to put on the graph
legendHIDAv = cell(beaconMax);  % To store the HIDs with an AV suffix
legendHIDSD = cell(beaconMax);  % To store the HIDs with an SD suffix

runningAv = zeros(maxRows, beaconMax);   %To store the running average values
runningSD = zeros(maxRows, beaconMax);  % To store the running SD values
badBeacons = zeros(beaconMax);

maxRunningAv(beaconMax) = -90; %initial maximum running average
minRunningAv(beaconMax) = -30; % initial minimum running average

maxRunningSD(beaconMax) = 0;   % initial maximum running SD
minRunningSD(beaconMax) = 100;  % initial minimum running SD

SDChangeThreshold = 1;
badSDThresh = 3;
%*****************************************

%**********Live Figure Setup and settings***************
% In order to have a live feed we need to assign plot lines to a variable
% or array of variables.  The lines can be adjusted using the set command
% and the data can be adjusted in a similar way (see below for this).
% displayNumber = 200;  % Number of points we wish to display
% tAxis = 1:displayNumber;    % vector for number of points to display
% measurePerBeacon = 3;   % The expected number of measurements per beacon
% axisMatrix = zeros(displayNumber,measurePerBeacon);        % Matrix for all of the results streams.  
% hold off
% rawRSSIplot = plot(tAxis,axisMatrix(:,1)); % This should be the line for the beacon raw RSSI
% hold on
% runningAVplot = plot(tAxis,axisMatrix(:,2));   % The plot for running average RSSI
% runningSDplot = plot(tAxis,axisMatrix(:,3));    % The plot for runnin SD
% 
% % Line spec settings for these plots
% set(rawRSSIplot, 'Color', 'b');
% set(runningAVplot, 'Color', 'r');
% set(runningSDplot, 'Color', 'k');
% ylim([-100,5]);

%******************************************

% %**********Get test file path*************
% suggestedInputPath = 'C:\Users\Mathew\Documents\MATLAB\BeaconTestData\*.csv';
% [inputFileName,inputFilePath] = uigetfile(suggestedInputPath);
% inputFullPath = strcat(inputFilePath, inputFileName);
% 
% if(strcmp(inputFileName,''))
%     %invalid file name, exit function
%     disp('No valid file.  Exiting')
% end
% %****************************************

%***************Static test file path*********
inputFilePath = 'C:\Documents and Settings\Administrator\My Documents\MATLAB\Beacon testing\BeaconTestData\';
inputFileName = 'beaconTester_testDataSet.csv';
inputFullPath = strcat(inputFilePath,inputFileName)
%***********************************************

%**************ONLY USED IN LIVE TRACKER
% % % Open comm port with correct settings
s = serial(strcat('COM1',serialPortNumber), 'baudrate', serialBaud);
s.InputBufferSize = serialInputBufferSize;
s.FlowControl = serialFlowControl;
s.Status;   % Diagnostic value to check port closed before opening
% % fopen(s);   %open the serial port
% % s.Status
% % 
% % if(~strcmp(s.Status,'open'))
% %    % Serial port failed to open.  Abort function
% %    disp('Port did not open')
% %    fclose(s);
% %    return
% % end
%*****************************************

%***********Read/Write File open**********
% Open input file for reading
% inputID = fopen(inputFullPath,'r'); %open the file that contains the test data
% if(inputID ==-1)
%     %failed to open file.  Display error and exit
%     disp('Input file failed to open.  exiting')
%     fclose(s);
%     return
% end

% Open logFile for writing (new File at each execution)
% logID = fopen(strcat(logFilePath,logFileName),'w');
% 
% if(logID == -1)
%     % File did not open.  Close everything and exit
%     disp('File did not open.  Closing serial and returning')
%     fclose(s);
%     fclose(inputID);
%     return
% end
%  %***************************************

 
% *************Get data from file*********************
rawData = csvread(inputFullPath);   % import all data from the test file
[rawRowCount, rawColCount] = size(rawData); % get the size of the data


% Parse data from file.  Note that this can be done in one hit when reading
% from a file.  It will need to be done on the fly when dragging from a
% stream.
for r=1:rawRowCount
    for c = 1:2:rawColCount-1   %This might need to be +/-1 depending
        % This inner loop steps in twos since the data is HID,RSSI format.
        % Modifiers to c must be applied to ensure that the correct field
        % is accessed.  This should be c = HID and c+1=RSSI
        
        % only need to collect HID's once
        if r == 1
            beaconHID(ceil(c/2)) = rawData(r,c);
        end     % End r==1 if
        
        % Regardless of column number, record the RSSI for that unit
        beaconRSSI(r,ceil(c/2)) = rawData(r,c+1);
    end     % End c rawColCount loop
end %End r rawRowCount for loop

% Write rawData to a logfile.  Note that this will need to be done
% line-by-line when getting data from serial.  When doing serial, recommend
% writting the whole line prior to parsing.
% for r = 1:rawRowCount
%     for c = 1:rawColCount
%         if c == rawColCount
%             %end of a column, prind the value with a linefeed
%             fprintf(logID,'%d\n',rawData(r,c));
%         else
%             %Any other part of the line.  do value followed by a comma
%             fprintf(logID,'%d,',rawData(r,c));
%         end     %end c=rawColCount if
%     end     % end c=1:rawColCount for
% end     %end r = 1:rawRowCount fo


% find out how many beacons are actually in the data
capturedBeacons = length(beaconHID);
% extract the names of the beacons for use in graphing
for b = 1:capturedBeacons
    legendHID(b) = {num2str(beaconHID(b))}; 
    legendHIDAv(b) = {strcat(num2str(beaconHID(b)),' AV')};
    legendHIDSD(b) = {strcat(num2str(beaconHID(b)),' SD')};
end
%****************End File read/serial extract***********

% ***************Process the data*****************

% Any recored zero value should be changed to the average of the last
% window's worth of data.  If the window is not filled then average over
% the last values that we have.  
for r = 1:rawRowCount
    for b = 1:capturedBeacons
        if beaconRSSI(r,b) == 0 && (r < runningAvWindow)
            % Haven't reached the running average window size.  Do the
            % average based on what we have
            beaconRSSI(r,b) = mean(beaconRSSI(1:r,b));
%             beaconRSSI(r,b) = -90;
        elseif beaconRSSI(r,b) == 0 
            % Row count greater than window, can apply running average to
            % the raw data.  NOTE:  make sure that this process does not
            % damage the assessment of the beacon.  Maybe record the number
            % of averages made to a beacon
        end
    end
end

for r = runningAvWindow+1:rawRowCount
    if r == 10
        r;
    end
    % Starting at the end of the first running average window, take the
    % mean of the window.  Potentially will need to do SD over this window
    % too.
    for b = 1:capturedBeacons
        % Do the running average and other stats for each beacon
        dataset=beaconRSSI(r-runningAvWindow:r-1,b);
        for d = 1:runningAvWindow
            if dataset(d) == 0
                if r == 6
                    dataset(d) = -65;
                else
                    dataset(d) = runningAv(r-runningAvWindow-1);
                end
            end
        end
        runningAv(r-runningAvWindow,b) = mean(dataset);        
        runningSD(r-runningAvWindow,b) = std(dataset);
        
        if maxRunningAv(b) < runningAv(r-runningAvWindow,b)&& (runningSD(r-runningAvWindow,b) < SDChangeThreshold)
             maxRunningAv(b) = runningAv(r-runningAvWindow,b);
        end
        if minRunningAv(b) > runningAv(r-runningAvWindow,b)&& (runningSD(r-runningAvWindow,b) < SDChangeThreshold)
             minRunningAv(b) = runningAv(r-runningAvWindow,b);
        end
        if( maxRunningSD(b) < runningSD(r-runningAvWindow,b)) 
             maxRunningSD(b) = runningSD(r-runningAvWindow,b);
        end
        if minRunningSD(b) > runningSD(r-runningAvWindow,b) 
             minRunningSD(b) = runningSD(r-runningAvWindow,b);
        end
        
%         newData = [beaconRSSI(r,b),runningAv(r), runningSD(r)]; %assemble the new data into a vector
%         axisMatrix(1,:) = [];   % Remove last row from graphed data
%         axisMatrix(displayNumber,:) = newData;  % Add new data to the end of the graph
%         set(rawRSSIplot,'YData',axisMatrix(:,1));   % Update raw plot
%         set(runningAVplot, 'YData', axisMatrix(:,2));   % Update running average plot
%         set(runningSDplot, 'YData', axisMatrix(:,3));   % Update running SD plot
        
%         drawnow;
    end % beacons for loop
end     % end running average rows for loop

% SD of the running av values
totalSD = std(runningAv);
deltaAv = maxRunningAv - minRunningAv;
deltaSD = maxRunningSD - minRunningSD

% Determine pass/fail and print to screen
fprintf(1,'Bad Beacons:\n');
for b = 1:capturedBeacons
    if deltaAv(b) == -60
        deltaAv(b) = 0
    else
        deltaAv(b)
    end
   
    if deltaSD(b) > badSDThresh
        % if the SD is outside the threshold then mark it as such on the
        % plot
        badBeacon(b) = 1;
        fprintf(1,'%d\n',beaconHID(b));
    end
end



% Graph the data %live
%plot(1:rawRowCount-5,runningAv,1:rawRowCount-5,runningSD,1:rawRowCount,beaconRSSI(:,1))
figure(1)
RSSIPlot = plot(1:rawRowCount,beaconRSSI(1:rawRowCount,:),'b-');
hold on
AVPlot = plot(1:rawRowCount,runningAv(1:rawRowCount,:),'r-');
SDPlot = plot(1:rawRowCount,runningSD(1:rawRowCount,:), 'k-');
hold off
legend([legendHID,legendHIDAv,legendHIDSD])


% When finished, close the ports and the file
% Close ports
fclose(s);  % Close the serial port
%s.Status
% if(exist('inputID') == 0)
%     fclose(inputID);
% end
% if(exist('logID') == 0)
%    fclose(logID);
% end
