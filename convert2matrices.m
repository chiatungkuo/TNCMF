%% This script processes the trip record matrix constructed from 'preprocess' and 
%% constructs 3 associated matrices as illustrated in the paper
%%  S: pickup location by dropoff location
%%  P: pickup location by pickup time
%%  D: dropoff location by dropoff time

% load('../trips.mat');       % load results from 'preprocess'

[PICKUP_LAT_COL, PICKUP_LONG_COL, PICKUP_TIME_COL, DROPOFF_LAT_COL, ...
    DROPOFF_LONG_COL, DROPOFF_TIME_COL] = deal(1, 2, 3, 4, 5, 6);

% set spatial and temporal granularities
LAT_LONG_GRANULARITY = 0.01;    %  0.01 degree for both latitude and longitude
TIME_GRANULARITY = 60;          % i.e. 1 minute
TIME_SPAN = 10 * 24 * 60;       % 10 days in minutes

% coarsen lat/long
data(:, PICKUP_LAT_COL) = round(data(:, PICKUP_LAT_COL)/LAT_LONG_GRANULARITY);
data(:, PICKUP_LONG_COL) = round(data(:, PICKUP_LONG_COL)/LAT_LONG_GRANULARITY);
data(:, DROPOFF_LAT_COL) = round(data(:, DROPOFF_LAT_COL)/LAT_LONG_GRANULARITY);
data(:, DROPOFF_LONG_COL) = round(data(:, DROPOFF_LONG_COL)/LAT_LONG_GRANULARITY);

% coarsen time
data(:, PICKUP_TIME_COL) = round(data(:, PICKUP_TIME_COL)/TIME_GRANULARITY);
data(:, DROPOFF_TIME_COL) = round(data(:, DROPOFF_TIME_COL)/TIME_GRANULARITY);

% throw away erroneous data points
% ignore trips spanning less than 1 minute or more than 3 hours (i.e. 180 minutes)
% ignore trips spanning more than 0.5 degrees of latitude & longitude
tripTime = data(:, DROPOFF_TIME_COL) - data(:, PICKUP_TIME_COL);
data = data(tripTime > 0 & tripTime <= 180, :);
latDist = abs(data(:, DROPOFF_LAT_COL) - data(:, PICKUP_LAT_COL));
data = data(latDist <= 0.5/LAT_LONG_GRANULARITY, :);
longDist = abs(data(:, DROPOFF_LONG_COL) - data(:, PICKUP_LONG_COL));
data = data(longDist <= 0.5/LAT_LONG_GRANULARITY, :);

% base of indices when put into matrices (i.e. first index corresponds to
% base + 1 in real data measurement)
baseLat = min([data(:, PICKUP_LAT_COL); data(:, DROPOFF_LAT_COL)]) - 1;
baseLong = min([data(:, PICKUP_LONG_COL); data(:, DROPOFF_LONG_COL)]) - 1;
baseTime = min([data(:, PICKUP_TIME_COL); data(:, DROPOFF_TIME_COL)]) - 1;

% determine each dimension size
dimLat = max([data(:, PICKUP_LAT_COL); data(:, DROPOFF_LAT_COL)]) - baseLat;
dimLong = max([data(:, PICKUP_LONG_COL); data(:, DROPOFF_LONG_COL)]) - baseLong;
dimTime = max([data(:, PICKUP_TIME_COL); data(:, DROPOFF_TIME_COL)]) - baseTime;

% record location data in index formats for construction of sparse matrices
P_LOC_INDEX = zeros(size(data, 1), 1);
D_LOC_INDEX = zeros(size(data, 1), 1);
P_TIME_INDEX = zeros(size(data, 1), 1);
D_TIME_INDEX = zeros(size(data, 1), 1);
VALUE = ones(size(data, 1), 1);

% dataIdxFormat = zeros(size(data, 1), 4);

for i = 1:size(data, 1)
    platIndex = data(i, PICKUP_LAT_COL);
    plongIndex = data(i, PICKUP_LONG_COL);
    
    dlatIndex = data(i, DROPOFF_LAT_COL);
    dlongIndex = data(i, DROPOFF_LONG_COL);
    
    ptimeIndex = data(i, PICKUP_TIME_COL);
    dtimeIndex = data(i, DROPOFF_TIME_COL);
    
    pIndex = latlong2ind(platIndex, plongIndex, dimLat, baseLat, baseLong);
    dIndex = latlong2ind(dlatIndex, dlongIndex, dimLat, baseLat, baseLong);
    
%     dataIdxFormat(i, :) = [pIndex, dIndex, ptimeIndex-baseTime, dtimeIndex-baseTime];
    
    P_LOC_INDEX(i) = pIndex;
    D_LOC_INDEX(i) = dIndex;
    
    P_TIME_INDEX(i) = ptimeIndex - baseTime;
    D_TIME_INDEX(i) = dtimeIndex - baseTime;
end
S = sparse(P_LOC_INDEX, D_LOC_INDEX, VALUE, dimLat*dimLong, dimLat*dimLong, size(data, 1));
P = sparse(P_LOC_INDEX, P_TIME_INDEX, VALUE, dimLat*dimLong, dimTime, size(data, 1));
D = sparse(D_LOC_INDEX, D_TIME_INDEX, VALUE, dimLat*dimLong, dimTime, size(data, 1));
