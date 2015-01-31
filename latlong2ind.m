function ind = latlong2ind(latInd, longInd, dimLat, baseLat, baseLong)
% Converts a given pair of latitude and longitude indices to an index in
% vectorized representation

ind =  (latInd - baseLat) + dimLat .* (longInd - baseLong - 1);
end