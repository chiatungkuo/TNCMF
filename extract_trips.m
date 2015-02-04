function trips = extract_trips(filepath, timeBreak)
% Extracts trip records from a given file containing GPS traces
% Inputs
%   filepath: string containing the path to the file
%   timeBreak: time in seconds such that two consecutive records further
%       apart are considered two separate trips
% Outputs
%   trips: n by 6 matrix where each row corresponds to a trip

SECONDS_PER_MINUTE = 60;
TIME_COLUMN = 4;

if ~exist('timeBreak', 'var'), timeBreak = 10 * SECONDS_PER_MINUTE; end

%% compute the time intervals between each pair of update records
data = load(filepath);
numRecords = size(data, 1);
intervals = zeros(numRecords-1, 1);
for row = 2:numRecords
    intervals(row-1) = data(row-1, TIME_COLUMN) - data(row, TIME_COLUMN);
end

%% break into smaller independent data sets if two updates are far apart
indicesToBreak = find(intervals > timeBreak);
trips = zeros(numRecords, 6);
numTrips = 0;

for i = 1:length(indicesToBreak)
    if i == 1
        prevIndex = 0;
    else
        prevIndex = indicesToBreak(i-1);
    end
    currIndex = indicesToBreak(i);
    subdata = data(prevIndex+1:currIndex, :);
    result = process(subdata);
    trips(numTrips+1:numTrips+size(result, 1), :) = result;
    numTrips = numTrips + size(result, 1);
end

% process the last sub-data
subdata = data(currIndex+1:end, :);
result = process(subdata);
trips(numTrips+1:numTrips+size(result, 1), :) = result;
numTrips = numTrips + size(result, 1);

trips = trips(1:numTrips, :);
end



function trips = process(data)
% subroutine that extracts trip records assuming the data set does not have
% discontinuous updates

LAT_COL = 1;
LONG_COL = 2;
FARE_COL = 3;
TIME_COL = 4;
numRecords = size(data, 1);
trips = zeros(numRecords, 6);
numTrips = 0;

% skip the initial 1's
curIndex = 1;
while data(curIndex, FARE_COL) == 1 && curIndex < numRecords, curIndex = curIndex + 1; end

% extract trip records
% notice that we record 'dropoff' before 'pickup' because the records are
% stored reverse chronologically (i.e. most recent on top of file)
curFare = 0;
while curIndex < numRecords
    fare = data(curIndex, FARE_COL);
    if fare == curFare
        curIndex = curIndex + 1;
        continue;
    elseif fare == 1 && curFare == 0
        dropoffUpdate = data(curIndex, [LAT_COL, LONG_COL, TIME_COL]);
        curFare = 1;
        curIndex = curIndex + 1;
        continue;
    else
        pickupUpdate = data(curIndex, [LAT_COL, LONG_COL, TIME_COL]);
        curFare = 0;
        trips(numTrips+1, :) = [pickupUpdate, dropoffUpdate];
        numTrips = numTrips + 1;
        curIndex = curIndex + 1;
    end
end
trips = trips(1:numTrips, :);
end
