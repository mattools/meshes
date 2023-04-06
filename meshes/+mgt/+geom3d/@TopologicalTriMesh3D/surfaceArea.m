function area = surfaceArea(obj)
% Surface area of this mesh, obtained by summing face areas.
%
% See Also
%   volume

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-27,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

% compute two direction vectors of each trinagular face, using the
% first vertex of each face as origin
v1 = obj.Vertices(obj.Faces(:, 2), :) - obj.Vertices(obj.Faces(:, 1), :);
v2 = obj.Vertices(obj.Faces(:, 3), :) - obj.Vertices(obj.Faces(:, 1), :);

% area of each triangle is half the cross product norm
% see also crossProduct3d in MatGeom
vn = zeros(size(v1));
vn(:) = bsxfun(@times, v1(:,[2 3 1],:), v2(:,[3 1 2],:)) - ...
    bsxfun(@times, v2(:,[2 3 1],:), v1(:,[3 1 2],:));
vn = sqrt(sum(vn .* vn, 2));

% sum up and normalize
area = sum(vn) / 2;