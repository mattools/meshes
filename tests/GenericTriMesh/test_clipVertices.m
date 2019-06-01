function test_suite = test_clipVertices(varargin)
%TEST_CLIPVERTICES  Test case for the file clipVertices
%
%   Test case for the file clipVertices

%   Example
%   test_clipVertices
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
% Test call of function without argument

oct = GenericTriMesh.createOctahedron();
box = [-0.5 2 -0.5 2 -0.5 2];
octC = clipVertices(oct, box);

testCase.assertEqual(3, vertexNumber(octC));
testCase.assertEqual(1, faceNumber(octC));


