function tests = test_addFace
% Test suite for the file addFace.
%
%   Test suite for the file addFace
%
%   Example
%   test_addFace
%
%   See also
%     addFace

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);


function test_AddSingleFace(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = mgt.geom3d.TopologicalTriMesh3D();

addVertex(mesh, [0 0 0]);
addVertex(mesh, [10 0 0]);
addVertex(mesh, [0 10 0]);

addFace(mesh, [1 3 2]);

assertEqual(testCase, faceCount(mesh), 1);
assertEqual(testCase, edgeCount(mesh), 3);


function test_Tetrahedron(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = mgt.geom3d.TopologicalTriMesh3D();

addVertex(mesh, [0 0 0]);
addVertex(mesh, [10 0 0]);
addVertex(mesh, [0 10 0]);
addVertex(mesh, [0 0 10]);

addFace(mesh, [1 3 2]);
addFace(mesh, [1 4 3]);
addFace(mesh, [1 2 4]);
addFace(mesh, [2 3 4]);

assertEqual(testCase, faceCount(mesh), 4);
assertEqual(testCase, edgeCount(mesh), 6);
