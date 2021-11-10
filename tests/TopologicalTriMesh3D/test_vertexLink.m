function tests = test_vertexLink
% Test suite for the file vertexLink.
%
%   Test suite for the file vertexLink
%
%   Example
%   test_vertexLink
%
%   See also
%     vertexLink

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

vertices = [...
    30 30 10; 50 30 10; 70 30 10; ...
    20 50 10; 40 50 10; 60 50 10; 80 50 10; ...
    30 70 10; 50 70 10; 70 70 10];
faces = [...
    1 2 5; 1 5 4; 2 3 6; 2 6 5; 3 7 6; ...
    4 5 8; 5 9 8; 5 6 9; 6 10 9; 6 7 10]; 
mesh = TopologicalTriMesh3D(vertices, faces);

link = vertexLink(mesh, 5);

assertEqual(testCase, vertexCount(link), 6);
assertEqual(testCase, edgeCount(link), 6);
assertEqual(testCase, faceCount(link), 0);

