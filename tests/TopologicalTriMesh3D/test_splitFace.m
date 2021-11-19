function tests = test_splitFace
% Test suite for the file splitFace.
%
%   Test suite for the file splitFace
%
%   Example
%   test_splitFace
%
%   See also
%     splitFace

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-19,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = TopologicalTriMesh3D();
v1 = addVertex(mesh, [10 10 10]);
v2 = addVertex(mesh, [30 10 10]);
v3 = addVertex(mesh, [20 40 10]);
f1 = addFace(mesh, [v1 v2 v3]);

splitFace(mesh, f1);

assertEqual(testCase, 4, vertexCount(mesh));
assertEqual(testCase, 6, edgeCount(mesh));
assertEqual(testCase, 3, faceCount(mesh));

