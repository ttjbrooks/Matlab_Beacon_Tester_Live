function [ tdNow ] = timeDateNow( )
% timeDateNow returns the date and time with no invalid chars so the name
% can be used to open files
%   [tdNow] = timeDateNow()
%   
%   Outputs:
%   tdNow = string containing the current date and time with no file
%   invalid characters
%
%   Inputs:
%   NONE

dateString = datestr(now);
tdNow = strrep(dateString,':','_');


end

