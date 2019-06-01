function test_suite = test_GenericTriMesh(varargin)
%TEST_GENERICTRIMESH  Test case for the file GenericTriMesh
%
%   Test case for the file GenericTriMesh

%   Example
%   test_GenericTriMesh
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-06-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% TEst constructor using sequential addition of elements

mesh = GenericTriMesh();
v1 = addVertex(mesh, [0 0 0]);
v2 = addVertex(mesh, [1 0 0]);
v3 = addVertex(mesh, [0 1 0]);
v4 = addVertex(mesh, [0 0 1]);
addFace(mesh, [v1 v3 v2]);
addFace(mesh, [v1 v2 v4]);
addFace(mesh, [v2 v3 v4]);
addFace(mesh, [v3 v1 v4]);

testCase.assertEqual(4, vertexNumber(mesh));
testCase.assertEqual(4, faceNumber(mesh));
computeEdges(mesh);
testCase.assertEqual(6, edgeNumber(mesh));
