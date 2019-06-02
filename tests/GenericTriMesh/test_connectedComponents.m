function test_suite = test_connectedComponents(varargin)
%TEST_VERTEXLINK  Test case for the file vertexLink
%
%   Test case for the file vertexLink

%   Example
%   test_vertexLink
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-06-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_TwoTetrahedra(testCase) %#ok<*DEFNU>
% Test call of function without argument

mesh = GenericTriMesh();
v1 = addVertex(mesh, [0 0 0]);
v2 = addVertex(mesh, [1 0 1]);
v3 = addVertex(mesh, [0 1 1]);
v4 = addVertex(mesh, [-0.7 -0.7 1]);
v5 = addVertex(mesh, [-1 0 -1]);
v6 = addVertex(mesh, [0 -1 -1]);
v7 = addVertex(mesh, [0.7 0.7 -1]);

addFace(mesh, [v1 v2 v3]);
addFace(mesh, [v1 v3 v4]);
addFace(mesh, [v1 v4 v2]);
addFace(mesh, [v2 v3 v4]);
addFace(mesh, [v1 v5 v6]);
addFace(mesh, [v1 v6 v7]);
addFace(mesh, [v1 v7 v5]);
addFace(mesh, [v5 v7 v6]);

ccList = connectedComponents(mesh);
testCase.assertEqual(2, length(ccList));
