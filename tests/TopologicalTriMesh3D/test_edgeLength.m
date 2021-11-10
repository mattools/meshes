function tests = test_edgeLength
% Test suite for the file edgeLength.
%
%   Test suite for the file edgeLength
%
%   Example
%   test_edgeLength
%
%   See also
%     edgeLength

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_AllEdges(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = TopologicalTriMesh3D();
addVertex(mesh, [0 0 0]);
addVertex(mesh, [10 0 0]);
addVertex(mesh, [0 10 0]);
addEdge(mesh, [1 2]);
addEdge(mesh, [1 3]);
addEdge(mesh, [2 3]);

els = edgeLength(mesh);

assertEqual(testCase, length(els), 3);
assertEqual(testCase, els, [10 10 10*sqrt(2)]');


function test_ChooseEdges(testCase) %#ok<*DEFNU>
% Test call of function without argument.

mesh = TopologicalTriMesh3D();
addVertex(mesh, [0 0 0]);
addVertex(mesh, [10 0 0]);
addVertex(mesh, [0 10 0]);
ie1 = addEdge(mesh, [1 2]);
ie2 = addEdge(mesh, [1 3]);
addEdge(mesh, [2 3]);

els = edgeLength(mesh, [ie1 ie2]);

assertEqual(testCase, length(els), 2);
assertEqual(testCase, els, [10 10]');



