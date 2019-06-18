function test_suite = test_flipEdge(varargin)
%TEST_FLIPEDGE  Test case for the file flipEdge
%
%   Test case for the file flipEdge

%   Example
%   test_flipEdge
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-06-18,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_TwoTriangles(testCase) %#ok<*DEFNU>
% boundary should contain four edges
vertices = [0 0 0; 0 4 0; 3 2 0;-3 2 0];
faces = [1 3 2;1 2 4];
mesh = GenericTriMesh(vertices, faces);

computeEdges(mesh);
ind = find(sum(ismember(mesh.Edges, [1 2]), 2) == 2); 

flipEdge(mesh, ind);

assertEqual(testCase, 2, sum(ismember(mesh.Edges(1,:), [3 4]), 2));
