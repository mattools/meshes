function test_suite = test_removeVertex(varargin)
%TEST_REMOVEVERTEX  Test case for the file removeVertex
%
%   Test case for the file removeVertex

%   Example
%   test_removeVertex
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-06-18,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument
obj = GenericTriMesh();
v1 = addVertex(obj, [5 4 3]);
v2 = addVertex(obj, [6 5 4]); %#ok<NASGU>
removeVertex(obj, v1);

testCase.assertEqual(1, vertexNumber(obj));

function test_WithFace(testCase) %#ok<*DEFNU>
% Test call of function without argument
vertices = [0 0 0; 0 4 0; 3 2 0;-3 2 0];
faces = [1 2 4];
mesh = GenericTriMesh(vertices, faces);
removeVertex(mesh, 3);
testCase.assertEqual(3, vertexNumber(mesh));
