function tests = test_collapseEdge
% Test suite for the file collapseEdge.
%
%   Test suite for the file collapseEdge
%
%   Example
%   test_collapseEdge
%
%   See also
%     collapseEdge

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
    1 2 5; 1 5 4; 2 3 6; 2 6 5; 3 7 6; ...
    4 5 8; 5 6 9; 5 9 8; 6 7 10; 6 10 9];
mesh = TopologicalTriMesh3D(vertices, faces);

% identify the edge (choose the one in the middle...)
ensureValidEdges(mesh);
indEdge = find(sum(ismember(mesh.Edges, [5 6]), 2) == 2);

collapseEdge(mesh, indEdge);

assertEqual(testCase, vertexCount(mesh), 9);
assertEqual(testCase, faceCount(mesh), 8);
