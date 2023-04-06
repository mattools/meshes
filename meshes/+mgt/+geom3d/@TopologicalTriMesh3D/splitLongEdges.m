function splitLongEdges(obj, maxLength)
%SPLITLONGEDGES Split all edges with length greater than given value.
%
%   splitLongEdges(MESH, MAXLENGTH)
%
%   Example
%   splitLongEdges
%
%   See also
%     splitEdge, collapseEdge, collapseSmallEdges
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2023-04-06,    using Matlab 9.13.0.2049777 (R2022b)
% Copyright 2023 INRAE.

% identifies edges with enough length 
[edgeLengths, inds] = edgeLength(obj);
inds = inds(edgeLengths > maxLength);

for iEdge = 1:length(inds)
    ind = inds(iEdge);
    splitEdge(obj, ind);
end
