function tests = test_removeInvalidVertices
% Test suite for the file removeInvalidVertices.
%
%   Test suite for the file removeInvalidVertices
%
%   Example
%   test_removeInvalidVertices
%
%   See also
%     removeInvalidVertices

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

% first create a simple mesh with two vertices close to each others
vertices = [...
    30 30; 50 30; 70 30; ...
    20 50; 40 50; 60 50; 80 50; ...
    30 70; 50 70; 70 70];
faces = [...
    1 2 6; 1 6 4; 2 3 6; 3 7 6; ...
    4 6 8; 6 9 8; 6 7 10; 6 10 9];
% faces = [...
%     1 2 6; 1 6 4; 2 3 6; 2 6 6; 3 7 6; ...
%     4 6 8; 6 6 9; 6 9 8; 6 7 10; 6 10 9];
mesh = TopologicalTriMesh3D(vertices, faces);
removeVertex(mesh, 5);
removeInvalidVertices(mesh);

assertEqual(testCase, 9, vertexCount(mesh));
assertEqual(testCase, 9, size(mesh.Vertices, 1));
