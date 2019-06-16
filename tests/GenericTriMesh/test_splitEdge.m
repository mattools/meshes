function test_suite = test_splitEdge(varargin)
%TEST_VERTEXLINK  Test case for the file vertexLink
%
%   Test case for the file vertexLink

%   Example
%   test_splitEdge
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-06-16,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_Octahedron(testCase) %#ok<*DEFNU>
% Test call of function without argument

oct = GenericTriMesh.createOctahedron();
splitEdge(oct, 4);

siz = size(oct);
assertEqual(testCase, siz, [7 15 10]);
