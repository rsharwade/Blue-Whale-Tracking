function [Tpred,xGrid,yGrid,nxGrid,nyGrid,range] = createTpred(staVec,xLim,yLim,xyBorder,dx,dy,velocity)
% Creates a 3D Matrix of predicted travel times
%
% Usage
%   [Tpred,xGrid,yGrid,nxGrid,nyGrid,range] = createTpred(staVec,xLim,yLim,xyBorder,dx,dy,velocity)
%
% Inputs
%   staVec   - Scalar station location structure
%   xLim     - Limits in x (km) for grid [min max] (excluding border)
%   yLim     - Limits in y (km) for grid [min max] (excluding border)
%   xBorder  - width of border to add around grid (km)
%   dx, dy   - Grid spacing
%   velocity - Velocity (km/s) for travel time calculations
%
% Outputs
%   Tpred - 3-D matrix of predicted travel times with dimensions (nx, ny, nStations)
%   xGrid - Vector of grid x values
%   yGrid - Vector of grid y values
%   nxGrid, nyGrid - x and y dimensions of grid
%   range - 3-D matrix corresponding to Tpred with distances to station from each point in grid 
%           Distances and travel times assume source at sea surface


xmin = floor((xLim(1)-xyBorder)/dx)*dx;
xmax = ceil((xLim(2)+xyBorder)/dx)*dx;
xGrid = xmin:dx:xmax;
nxGrid = length(xGrid);

ymin = floor((yLim(1)-xyBorder)/dy)*dy;
ymax = ceil((yLim(2)+xyBorder)/dy)*dy;
yGrid = ymin:dy:ymax;
nyGrid = length(yGrid);

xMat = repmat(xGrid',1,nyGrid);
yMat = repmat(yGrid,nxGrid,1);

nsta = length(staVec.name);

Tpred = zeros(nxGrid,nyGrid,nsta);
if nargout>=5
  range = zeros(nxGrid,nyGrid,nsta);
end

for ista = 1:nsta
  Tpred(:,:,ista) = sqrt((xMat-staVec.x(ista)).^2 + ...
                         (yMat-staVec.y(ista)).^2 + ...
                         staVec.z(ista)^2)./velocity;
  if nargout>=5
    range(:,:,ista) = sqrt((xMat-staVec.x(ista)).^2 + ...
                           (yMat-staVec.y(ista)).^2 + ...
                           staVec.z(ista)^2);
  end
end

