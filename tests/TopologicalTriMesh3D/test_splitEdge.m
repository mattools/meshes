function tests = test_splitEdge
% Test suite for the file splitEdge.
%
%   Test suite for the file splitEdge
%
%   Example
%   test_splitEdge
%
%   See also
%     splitEdge

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
vf1 = addVertex(mesh, [20 30 10]);
vf2 = addVertex(mesh, [20 10 10]);
e0 = addEdge(mesh, v1, v2);
addFace(mesh, [v1 v2 vf1]);
addFace(mesh, [v1 vf2 v2]);

splitEdge(mesh, e0);

assertEqual(testCase, 5, vertexCount(mesh));
assertEqual(testCase, 8, edgeCount(mesh));
assertEqual(testCase, 4, faceCount(mesh));


