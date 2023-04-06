function h = drawFaceNormals(obj, varargin)
% Draw the normal vector of each mesh face.
%
%   output = truc(input)
%
%   Example
%   truc
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

% compute vector data
c = faceCentroids(obj);
n = faceNormals(obj);

% display an arrow for each normal
hq = quiver3(c(:,1), c(:,2), c(:,3), n(:,1), n(:,2), n(:,3));

% format output
if nargout > 0
    h = hq;
end
