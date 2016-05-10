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

%******************************************

% ******Define local variables***********
testCycle = 10;  % The number of packets to receive from the GL before test is complete.

beaconMax = 30; % The maximum number of beacons expected in a packet

% Define the max rows based on either the number of test cycles to run or
% by constant
if exist('testCycle', 'var')
    % testCycle exists, set the max number of rows the same
    maxRows = testCycle; % The maximum number of rows that we expect to receive (note that if the actual number is larger then the system will auto grow the array
else
    maxRows = 100;
end

beaconHID = zeros(1,beaconMax); % An array for the HID
beaconRSSI = zeros(maxRows,beaconMax); % An array for the RSSI (note, will grow in a loop) 
firstEmptyBeaconCell = 1;   %The index of the first empty cell in the HID array.  This number should only ever max out at beaconMax + 1 but this will not affect operation of this program

rowCounter = 1; % A counter for the number of rows that we have counted.  Always assume that there will be one row to deal with



serialPortNumber = '78';   % The number of the target comm port
serialBaud = 115200;    % Target baud rate for the serial
serialInputBufferSize = 4048;   % defined size for the input buffer
serialFlowControl = 'hardware';  % defined flow control type

serialDataLine = '';    % Variable to store a line of serial data from the port
serialCount = 0;    % Char count from the serial stream (should be constant ish)
serialErrMsg = 0;   % Storage for error message from the serial stream (may use this as a means of terminating program)

dataLine = cell(1,2*beaconMax+1);  % Array to hold the delimited serial data.  Pre allocate based on beaconMax
dataLineNum = zeros(1,2*beaconMax+1); % Array to hold the delimited data when converted to numbers.
dataLen = 0;    % length counter for the serial data.  Confirm if packet always at max beacon length

% modC = 0;   % Used to store modulus of column count.  Will determine odd or even columns
lineHID = zeros(1,beaconMax); % The HID values captured in the current data packet
lineRSSI = zeros(1,beaconMax);    % The RSSI values captured in the current data packet


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
legendHID = cell(1,beaconMax);    % To store the names of HIDs to put on the graph
% legendHIDAv = cell(1,beaconMax);  % To store the HIDs with an AV suffix
% legendHIDSD = cell(1,beaconMax);  % To store the HIDs with an SD suffix

legendHID(:) = {'0'};   % Initialize the legend HID cells to zero (the default for a missing beacon)

legendLocation = [0.2 1 3 16];

correctionCount = zeros(1,beaconMax);   % To store the number of times that a rawRSSI value needed to be corrected.
differentBeacons = false;

runningAv = zeros(maxRows, beaconMax);   %To store the running average values
runningSD = zeros(maxRows, beaconMax);  % To store the running SD values

totalAv = zeros(1,beaconMax);   % To store the total average RSSI for the entire test (per beacon)
totalSD = zeros(1,beaconMax);   % To store the total SD for the entire test (per beacon)

badBeacons = zeros(1,beaconMax);

maxRunningAv = zeros(1,beaconMax) -90; %initial maximum running average
minRunningAv = zeros(1,beaconMax) -30; % initial minimum running average

maxRunningSD = zeros(1,beaconMax);   % initial maximum running SD
minRunningSD = zeros(1,beaconMax) + 100;  % initial minimum running SD

%***********Thresholds for detecting bad beacons**************
SDChangeThreshold = 1;
badSDThresh = 3;
deltaSDThresh = 20;
%*************************************************************

%*****************Graphing options*****************************
drawDivider = 3;    % Figures redrawn every 'drawDivider' rows.
%*************************************************************

%**********Live Figure Setup and settings***************
% In order to have a live feed we need to assign plot lines to a variable
% or array of variables.  The lines can be adjusted using the set command
% and the data can be adjusted in a similar way (see below for this).

% Define the draw details for the figure.  These details are as follows:
% Row 1 (outer position base): [left bottom width height]
% note that it is intended that units be cm
figPos = zeros(2,4);
figPos(1,:) = [0 2 30 20];
figPos(2,:) = [figPos(1,1) + 0.5, figPos(1,2)+0.5, figPos(1,3)-1,figPos(1,4)];
displayNumber = maxRows;  % Number of points we wish to display
tAxis = 1:displayNumber;    % vector for number of points to display
measurePerBeacon = 3;   % The expected number of measurements per beacon
axisMatrix = zeros(displayNumber,beaconMax,measurePerBeacon);        % Matrix for all of the results streams. 
axisMatrix(:,:,1) = axisMatrix(:,:,1) - 90;   % Initialize the rawRSSI to minimum
axisMatrix(:,:,2) = axisMatrix(:,:,2) - 90;   % Initialize the runningAv to minimum


rawRSSIFig = figure(1);
hold off
subplot(3,1,1)
rawRSSIPlot = plot(tAxis,axisMatrix(:,:,1)); % This should be the line for the beacon raw RSSI
hold on
set(rawRSSIFig, 'Units', 'centimeters');
set(rawRSSIFig, 'OuterPosition', figPos(1,:));
title('Raw RSSI')
ylim([-100,50]);


subplot(3,1,2)
runningAvPlot = plot(tAxis,axisMatrix(:,:,2));   % The plot for running average RSSI
title('Running Average RSSI')



subplot(3,1,3)
runningSDPlot = plot(tAxis,axisMatrix(:,:,3));   % The plot for running average RSSI
title('Running SD')


% Line spec settings for these plots
% set(rawRSSIPlot, 'Color', 'b');
% set(runningAVplot, 'Color', 'r');
% set(runningSDplot, 'Color', 'k');
% ylim([-100,5]);

%******************End live figure settings************************

%***** Open comm port with correct settings**********************
s = serial(strcat('COM',serialPortNumber), 'baudrate', serialBaud);
s.InputBufferSize = serialInputBufferSize;
s.FlowControl = serialFlowControl;
%s.Status   % Diagnostic value to check port closed before opening
fopen(s);   %open the serial port
s.Status;

if(~strcmp(s.Status,'open'))
   % Serial port failed to open.  Abort function
   disp('Port did not open.  Closing all open connections and exiting')
   instrreset; % Close any open serial ports.
   return   % Abort function
end
%***************End com port config*************************

 
% % % ******* Open logFile for writing (new File at each execution)****
% % %logID = fopen(strcat(logFilePath, logFileName),'w');
% % logID = fopen(strcat(logFilePath,logFileName),'w');
% % logID
% % 
% % if(logID == -1)
% %     % File did not open.  Close everything and exit
% %     disp('File did not open.  Closing serial and returning')
% %     fclose(s);
% % %    fclose(logID);
% %     return
% % end

% log file open diagnostic
 %fprintf(logID,'Hello world'); Code working
 %******************************************************************
 


%  Loop through the required number of test cycles.  Note that pressing
%  'q' with the figure selected will end execution correctly.  Terminal
%  will display 'Exiting Beacon Tester' when exiting.
for cycleCount = 1:testCycle
 
   
    %******* Get packet from Serial stream ******
    [serialDataLine, serialCount, serialErrMsg] = fgetl(s);

    % % %Serial read diagnostic code
    % % serialDataLine
    % % serialCount
    % % serialErrMsg
    % % %End serial diagnostic code

    if serialErrMsg ~= ''
        %There was an error in reading the serial data line.  Print the error
        %and exit
        disp(serialErrMsg)
        disp('Closing all open ports')
        instrreset;
        return
    end

    % **********End get serial data*****************************
    
    %************Parse data*************************************
    % Parse the string into HID and RSSI (keep index the same for both)
    % Note strsplit was introduced in 2013a.  Running on previous versions will
    % require the standalone .m file.
    dataLine = strsplit(serialDataLine,',');
    dataLineNum = cell2mat(cellfun(@str2num,dataLine,'UniformOutput', false));

    % Separate data into HID and rawRSSI values.
    dataLen = length(dataLineNum);  %find the length of the packet
    
   % Use indexes steped in two to get the HID and RSSI from the packet.
   % Note that this array is variable length and we must account for this
   % when consolodating the values
    lineHID = dataLineNum(1:2:dataLen-1);   % Grab every other value and put it in line HID
    lineRSSI = dataLineNum(2:2:dataLen);    % Grab every other value (starting at 2) and put it into the lineRSSI
    

    % compare the lineHID to current HID values and update if necessary.
    % Note order of HIDs
    lineCapBeacons = length(find(lineHID ~= 0));    % Number of captured beacons in this packet
    differentBeacons = false;    
    
    for b = 1:lineCapBeacons
        % For each of the beacons captured in this line, find the location
        % in the beaconHID list then apply the RSSI value to the same
        % column in the beaconRSSI list at the current row.
        indexHID = find(beaconHID == lineHID(b));   % get the index of current HID in the current line
        newBeacon = isempty(indexHID);
        if newBeacon && (firstEmptyBeaconCell <= beaconMax)
            % The HID does not appear in our list.  Add it along with the
            % RSSI for the current row
            beaconHID(firstEmptyBeaconCell) = lineHID(b); % Transfer the HID to the next empty space
            differentBeacon = true; % the beacon list has been updated.
            beaconRSSI(rowCounter,firstEmptyBeaconCell) = lineRSSI(b); % Transfer the RSSI to the next empty space at the current row
            firstEmptyBeaconCell = firstEmptyBeaconCell + 1;    % Increment the next empty cell index
        elseif newBeacon && (firstEmptyBeaconCell > beaconMax)
            % The HID is not in our list but the list is full.  Display a
            % message and continue
            fprintf(1,'WARNING: More than %d beacons detected.  Increase beaconMax or clean up radio area.\n',beaconMax);
        else
            % The HID does appear in our list.  Add the RSSI at the correct
            % row and with the correct index
            beaconRSSI(rowCounter,indexHID) = lineRSSI(b);
        end
    end      
    
    %***********End Parsing*******************************************
    

    % **************Process the data**********************************
    
    % We now have raw data and corresponding HIDs.  We can calculate
    % running statistics.
    
    % Correct any incoming 0 RSSI values by taking the average over the
    % last window period.  If we are less than the window then use the
    % window size we have.  If there are zeros to the end of the window,
    % assume that the RSSI is -90
    
    beaconCount = length(find(beaconHID > 0));
    for b = 1:beaconCount
        % Check for a current 0 RSSI value
        if beaconRSSI(rowCounter,b) == 0 && rowCounter == 1
            % If we are at the first row then there are no previous values,
            % set RSSI to -90
            beaconRSSI(rowCounter,b) = -90;
            correctionCount(b) = correctionCount(b) + 1;
        elseif (beaconRSSI(rowCounter,b) == 0) && (rowCounter < runningAvWindow)
            % The value is zero, there are non-zero values before hand but
            % we have not reached enough rows to perform a full average.
            % Use the number of values that we do have
            beaconRSSI(rowCounter,b) = mean(beaconRSSI(1:rowCounter,b));    % the mean of the values that we do have
            correctionCount(b) = correctionCount(b) + 1;     % Increment the number of times we have had to correct the data
        elseif (beaconRSSI(rowCounter,b) == 0) && (rowCounter >= runningAvWindow)
            % RSSI zero and we have enough to perform an average.  Do a
            % running average estimation of the current RSSI
            beaconRSSI(rowCounter,b) = mean(beaconRSSI(rowCounter-runningAvWindow+1:rowCounter,b));
            correctionCount(b) = correctionCount(b) + 1;    
        else
            % Value must not be zero, do nothing            
        end
    end
    
    % Perform the running statistics
    if rowCounter < runningAvWindow +1
        % not enough data to do a full running stat, instead do with the
        % data we have.
        dataSet = beaconRSSI(1:rowCounter,:);
    else
        % There is enough data to perform the full window statistics.
        dataSet = beaconRSSI(rowCounter-runningAvWindow:rowCounter,:);
    end
    runningAv(rowCounter,:) = mean(dataSet);
    runningSD(rowCounter,:) = std(dataSet);
    
    % Perform the total statistics
    totalAv = mean(beaconRSSI(1:rowCounter,:));
    totalSD = std(beaconRSSI(1:rowCounter,:));
    
    
   if (rowCounter > runningAvWindow)
       for b = 1:beaconCount
            % Need to find a way to do this without the loop.  Somnething
            % weird is going on here too.
            if maxRunningAv(b) < runningAv(rowCounter-runningAvWindow,b)%&& (runningSD(rowCounter-runningAvWindow,b) < SDChangeThreshold)
                 maxRunningAv(b) = runningAv(rowCounter-runningAvWindow,b);
            end
            if minRunningAv(b) > runningAv(rowCounter-runningAvWindow,b)%&& (runningSD(rowCounter-runningAvWindow,b) < SDChangeThreshold)
                 minRunningAv(b) = runningAv(rowCounter-runningAvWindow,b);
            end
            if( maxRunningSD(b) < runningSD(rowCounter-runningAvWindow,b)) 
                 maxRunningSD(b) = runningSD(rowCounter-runningAvWindow,b);
            end
            if minRunningSD(b) > runningSD(rowCounter-runningAvWindow,b) 
                 minRunningSD(b) = runningSD(rowCounter-runningAvWindow,b);
            end
        end
   end
    
    
    %*************End data processing*********************************

    % ************* Graph the data live*******************************
    % This section takes the statistics calculated in the previous section
    % and applies it to the displayed graph.  It will be necessary to
    % manage the graphs and displayed statistics so that it is not too
    % cumbersome.
    
    % Remove first line of all sections of the axisMatrix (i.e. purge the
    % front of the queue.
    axisMatrix(1,:,:) = []; % Remove the line by nulling
    
    % Add data to the end of queue from our measurements
    axisMatrix(displayNumber,:,1) = beaconRSSI(rowCounter,:);   % Add the rawRSSI for this row to the graphing matrix
    axisMatrix(displayNumber,:,2) = runningAv(rowCounter,:);    % Add the runningAv for this row to the graphing matrix
    axisMatrix(displayNumber,:,3) = runningSD(rowCounter,:);    % Add the runningSD for this row to the graphing matrix
    axisMatrix(displayNumber,:,4) = minRunningSD(:);    % Track the minimum Running SD.  Note, this might not work
    axisMatrix(displayNumber,:,5) = maxRunningSD(:);    % Tracl the maximum Running SD. Note, this might not work.
    
    % Update the 'YData' parts of the figure. may need to be done in a loop
    for b=1:beaconCount
        % update the yData for each of the beacon plots one beacon at a
        % time
        set(rawRSSIPlot(b), 'YData' , axisMatrix(:,b,1));  % Update the rawRSSI to the figure
        set(runningAvPlot(b), 'YData', axisMatrix(:,b,2)); % Update the running Av Plot to the figure
        set(runningSDPlot(b), 'YData', axisMatrix(:,b,3));     % Update the runningSD plot to the figure
        if differentBeacon == true
            % There is a different beacon in the list therefore we must
            % change the legend
            legendHID(b) = {num2str(beaconHID(b))}; %Update the legendHID cell array
        end
    end
    
    % Update the legend of the plots
    if differentBeacon == true
        % There are different beacons in the list, update the legend
%         figure(rawRSSIFig)
        subplot(3,1,3)
        sdLeg = legend(legendHID, 'Location', 'northeastoutside');
        set(sdLeg, 'Units', 'centimeters')
        set(sdLeg, 'Position', legendLocation)

%         if rowCounter == 1
%             subplot(3,1,1)
%             avLeg = legend(legendHID, 'Location', 'northeastoutside');
%             set(avLeg,'visible','off');
% 
%             subplot(3,1,2)
%             rawLeg = legend(legendHID, 'Location', 'northeastoutside');
%             set(rawLeg, 'visible', 'off');
%         end

    end
    
    % redraw the figure
    if mod(rowCounter,drawDivider) == 0
        drawnow;
    end
    
    %************ End Live Graph *************************************
    

    %*** Determine pass/fail IDs and display data to user*************
    % This section needs Lucia's input to determine how we assess pass/fail
    % for beaons.  
    % For the moment, assume that if the SD is above a given
    % level that it is a bad beacon.  Need to determine a way to switch
    % badness on and off to account for bad readings at the beginning of
    % measurements
    
    % Only check for bad beacons after measuremetns have stabilized (this
    % value may need alteration
    if rowCounter > 5
        deltaSD = maxRunningSD - minRunningSD;
        badBeacons(deltaSD > deltaSDThresh) = 1;
    %     % The following code might be redundant
    %     for b = 1:beaconCount
    %         % Do all of the checks for a given beacon to decide if it's bad
    %         if deltaSD(b) > deltaSDThresh
    %             % This beacon is bad, mark it as such
    %             badBeacons(b) = 1;   
    %         end
    %     end
    %       % End potentially redundant code
    end

    % Display the bad beacons to the screen and change the legend on the SD
    % plot
    badLoc = find(badBeacons == 1);
    if ~isempty(badLoc)
        % This code should execute only if there are bad beacons
        fprintf(1,'Bad Beacons:\n');
        disp(beaconHID(badBeacons == 1));
        set(runningSDPlot(badLoc),'Marker','x');
    end
    %******************************************************************

%     Decide if the user has pressed the escape command
    figure(rawRSSIFig)
    figChar = double(get(rawRSSIFig,'CurrentChar'));
    if(figChar == 113)
        % quit char has been given.  Break the loop.  Else continue
        break;
    end
    
    % Done with this row of the stream.  Increment the row counter
%     fprintf(1,'Row number = %d\n',rowCounter); % Diagnostic progress counter
    rowCounter = rowCounter + 1;

    % end
end % end test for

drawnow

% When finished, close the ports and the file
% Close ports
fclose(s);  % Close the serial port
s.Status;
% fclose(logID);
disp('Exiting Beacon Tester');
