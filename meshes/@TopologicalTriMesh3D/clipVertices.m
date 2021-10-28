function res = clipVertices(obj, box)
% Clip the mesh by retaining only vertices within the box.
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

% clip the vertices
% get bounding box limits
xmin = box(1); xmax = box(2);
ymin = box(3); ymax = box(4);
zmin = box(5); zmax = box(6);

% compute indices of points inside visible area
xOk = obj.Vertices(:,1) >= xmin & obj.Vertices(:,1) <= xmax;
yOk = obj.Vertices(:,2) >= ymin & obj.Vertices(:,2) <= ymax;
zOk = obj.Vertices(:,3) >= zmin & obj.Vertices(:,3) <= zmax;

% select vertices
inds = find(xOk & yOk & zOk);
newVertices = obj.Vertices(inds, :);

% create index array for face indices relabeling
refInds = zeros(1, length(xOk));
for i = 1:length(inds)
    refInds(inds(i)) = i;
end

% select the faces with all vertices within the box
indFaces = sum(~ismember(obj.Faces, inds), 2) == 0;
newFaces = refInds(obj.Faces(indFaces, :));

res = TopologicalTriMesh3D(newVertices, newFaces);
