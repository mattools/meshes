function tests = test_isBoundaryEdge
% Test suite for the file isBoundaryEdge.
%
%   Test suite for the file isBoundaryEdge
%
%   Example
%   test_isBoundaryEdge
%
%   See also
%     isBoundaryEdge

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = TopologicalTriMesh3D();
v1 = addVertex(mesh, [ 0 10 10]);
v2 = addVertex(mesh, [10  0 10]);
v3 = addVertex(mesh, [10 20 10]);
v4 = addVertex(mesh, [20 10 10]);
addFace(mesh, [v1 v2 v3]);
addFace(mesh, [v2 v4 v3]);

edge12 = findEdgeIndex(mesh, [v1 v2]);
edge23 = findEdgeIndex(mesh, [v2 v3]);
edge24 = findEdgeIndex(mesh, [v2 v4]);

assertTrue(testCase, isBoundaryEdge(mesh, edge12));
assertFalse(testCase, isBoundaryEdge(mesh, edge23));
assertTrue(testCase, isBoundaryEdge(mesh, edge24));
