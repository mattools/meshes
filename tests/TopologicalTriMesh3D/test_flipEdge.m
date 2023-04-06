function tests = test_flipEdge
% Test suite for the file flipEdge.
%
%   Test suite for the file flipEdge
%
%   Example
%   test_flipEdge
%
%   See also
%     flipEdge

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-19,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = mgt.geom3d.TopologicalTriMesh3D();
v1 = addVertex(mesh, [10 20 10]);
v2 = addVertex(mesh, [30 20 10]);
v3 = addVertex(mesh, [20 30 10]);
v4 = addVertex(mesh, [20 10 10]);
e0 = addEdge(mesh, v1, v2);
addFace(mesh, [v1 v2 v3]);
addFace(mesh, [v1 v4 v2]);

flipEdge(mesh, e0);

assertEqual(testCase, 4, vertexCount(mesh));
assertEqual(testCase, 5, edgeCount(mesh));
assertEqual(testCase, 2, faceCount(mesh));

assertEqual(testCase, 2, length(mesh.VertexEdges{v1}));
assertEqual(testCase, 2, length(mesh.VertexEdges{v2}));
assertEqual(testCase, 3, length(mesh.VertexEdges{v3}));
assertEqual(testCase, 3, length(mesh.VertexEdges{v4}));

