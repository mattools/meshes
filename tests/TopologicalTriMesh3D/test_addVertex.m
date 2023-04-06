function tests = test_addVertex
% Test suite for the file addVertex.
%
%   Test suite for the file addVertex
%
%   Example
%   test_addVertex
%
%   See also
%     addVertex

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = mgt.geom3d.TopologicalTriMesh3D();

addVertex(mesh, [0 0 0]);
addVertex(mesh, [10 0 0]);
addVertex(mesh, [0 10 0]);
addVertex(mesh, [0 0 10]);

assertEqual(testCase, vertexCount(mesh), 4);


