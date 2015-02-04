%% Extract trip records from each taxi and put all into a matrix
path = './';        % set this to the directory containing those taxi GPS traces
files = dir(path);

MAX_RECORDS = 1e7;  % some large number to initialize array size
MAX_CABS = 536;

% each record consists of 6 attributes (see README)
data = zeros(MAX_RECORDS, 6);
numRecords = 0;

for i = 1:length(files)
    if strfind(files(i).name, 'new') == 1   % this looks for all file names with 'new' which were the default naming from data source
        trips = extract_trips(strcat(path, files(i).name));
        numTrips = size(trips, 1);
        data(numRecords+1:numRecords+numTrips, 1:6) = trips;
        numRecords = numRecords + numTrips;
    end
end

data = data(1:numRecords, :);
