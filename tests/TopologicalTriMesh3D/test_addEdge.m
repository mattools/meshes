function tests = test_addEdge
% Test suite for the file addEdge.
%
%   Test suite for the file addEdge
%
%   Example
%   test_addEdge
%
%   See also
%     addEdge

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test adding valid edges

mesh = mgt.geom3d.TopologicalTriMesh3D();
addVertex(mesh, [0 0 0]);
addVertex(mesh, [10 0 0]);
addVertex(mesh, [0 10 0]);
addVertex(mesh, [0 0 10]);

addEdge(mesh, [1 2]);
addEdge(mesh, [2 3]);
addEdge(mesh, [3 1]);
addEdge(mesh, [1 4]);
addEdge(mesh, [2 4]);
addEdge(mesh, [3 4]);

assertEqual(testCase, edgeCount(mesh), 6);
