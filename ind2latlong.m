function locIndex = ind2latlong(index, dimLat, baseLat, baseLong)
% Converts a given index (in vectorized format) to a pair of latitude and
% longitude indices

latIndex = mod(index, dimLat);
if latIndex == 0, latIndex = dimLat; end
longIndex = floor(index / dimLat) + 1;
locIndex = [latIndex + baseLat, longIndex + baseLong];
end
