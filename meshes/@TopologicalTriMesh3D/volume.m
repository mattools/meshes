function vol = volume(obj)
% (signed) volume enclosed by this mesh.
%
% See Also
%   surfaceArea

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-27,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

% initialize an array of volume
nFaces = size(obj.Faces, 1);
vols = zeros(nFaces, 1);

% Shift all vertices to the mesh centroid
centroid = mean(obj.Vertices, 1);

% compute volume of each tetraedron
for iFace = 1:nFaces
    % consider the tetrahedron formed by face and mesh centroid
    tetra = obj.Vertices(obj.Faces(iFace, :), :);
    tetra = bsxfun(@minus, tetra, centroid);
    
    % volume of current tetrahedron
    vols(iFace) = det(tetra) / 6;
end

vol = sum(vols);