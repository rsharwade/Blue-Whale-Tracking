function  detVec = compileDetections(detections)
% Merges detections for 1 day into single element structure
%
% Usage
%   detVec = compileDetections(detections)
% 
% Inputs
%   detections - Vector structure of detections from Erik/Rose
%
% detVec
%   Single structure with vector fields of all detections sorted by time
%   Fields comprise stationName, peakTime, peakAmp, duration, strength = duration * peakAmp

detVec.stationName = [];
detVec.peakTime = [];
detVec.peakAmp = [];
detVec.duration = [];

for i = 1:length(detections)
  if ~isempty(detections(i).peakTime)
    n = length(detections(i).peakTime);
    if isfield(detections,'stationID')
      detVec.stationName = [detVec.stationName; repmat({detections(i).stationID},n,1)];
    else
      detVec.stationName = [detVec.stationName; repmat({detections(i).stationName},n,1)];
    end
    detVec.peakTime = [detVec.peakTime; detections(i).peakTime];
    detVec.peakAmp = [detVec.peakAmp; detections(i).peakAmp];
    detVec.duration = [detVec.duration; 86400*(detections(i).detEnd-detections(i).detStart)];
  end
end

[detVec.peakTime,i] = sort(detVec.peakTime);
detVec.stationName = detVec.stationName(i);
detVec.peakAmp = detVec.peakAmp(i);
detVec.duration = detVec.duration(i);
detVec.strength = detVec.duration.*detVec.peakAmp;
