function res = smooth(obj, varargin)
% Apply smoothing on mesh.
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
% Created: 2021-10-27,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.
    

% determine number of iterations
nIter = 1;
if ~isempty(varargin)
    nIter = varargin{1};
end

% compute adjacency matrix,
% result is a Nv-by-Nv matrix with zeros on the diagonal
adj = adjacencyMatrix(obj);

% ensure the size of the matrix is Nv-by-Nv
% (this can not be the case if some vertices are not referenced)
nv = vertexCount(obj);
if size(adj, 1) < nv
    adj(nv, nv) = 0;
end

% Add "self adjacencies"
adj = adj + speye(nv);

% weight each vertex by the number of its neighbors
w = spdiags(full(sum(adj, 2).^(-1)), 0, nv, nv);
adj = w * adj;

% do averaging to smooth the field
v2 = obj.Vertices;
for k = 1:nIter
    v2 = adj * v2;
end

% create a new mesh with same faces
res = TopologicalTriMesh3D(v2, obj.Faces);
