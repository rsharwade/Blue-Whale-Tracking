function staVec = compileStations(station)
% Merges vector stations structure into a scalar structure with a matrix of station spacings
%
% Usage
%   staVec = compileStations(station);
%
% Input
%   station - Structure of station information.  One element per station with fields 
%             name, lat, lon and elev (m) fields
%               ** Station should track the on/off time for the station (when it recorded data)
%
% Output
%   staVec  - Scalar structure of station information with 
%             Vectors of name, lat, lon, elev (km), x (km), y (km)
%               *** Need to worry about what happens when stations are in multiple utm zones
%             Square matrix spacing with distance between stations in km 


for i=1:length(station)
  if length(station(i).stationName)==3;
    station(i).stationName = [station(i).stationName 'A'];
  end
end

staVec.name = {station.stationName}';
staVec.lat = [station.lat]';
staVec.lon = [station.long]';
staVec.elev = [station.elev];

[staVec.x,staVec.y,staVec.utmzone] = wgs2utm(staVec.lat,staVec.lon,9,'N');
staVec.x = staVec.x/1000;
staVec.y = staVec.y/1000;
staVec.z = staVec.elev/1000;

staVec.spacing = zeros(length(station));
for i = 1:length(station)
  staVec.spacing(:,i) = sqrt((staVec.x-staVec.x(i)).^2 + (staVec.y-staVec.y(i)).^2);
end



