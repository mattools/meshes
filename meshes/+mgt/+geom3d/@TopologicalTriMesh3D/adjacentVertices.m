function neighs = adjacentVertices(obj, vertexIndex)
% Retrieve adjacent vertex indices of a given vertex.
%
%   NEIGHS = adjacentVertices(MESH, VERT_IDX)
%
%   Example
%   adjacentVertices
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2022-07-19,    using Matlab 9.12.0.1884302 (R2022a)
% Copyright 2022 INRAE.


neighs = unique(obj.Faces(obj.VertexFaces{vertexIndex}, :));
neighs(neighs == vertexIndex) = [];
