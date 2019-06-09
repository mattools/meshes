function test_suite = test_isBoundaryVertex(varargin)
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
% Created: 2019-06-09,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_Otahedron(testCase) %#ok<*DEFNU>
% Test call of function without argument

oct = GenericTriMesh.createOctahedron;
removeFace(oct, 8);

bv = false(1, 6);
for i = 1:6
    bv(i) = isBoundaryVertex(oct, i);
end

testCase.assertEqual(3, sum(bv));
