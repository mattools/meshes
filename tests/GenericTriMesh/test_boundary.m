function test_suite = test_boundary(varargin)
%TEST_BOUNDARY  Test case for the file boundary
%
%   Test case for the file boundary

%   Example
%   test_boundary
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-06-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

test_suite = functiontests(localfunctions);

function test_TwoTriangles(testCase) %#ok<*DEFNU>
% boundary should contain four edges
vertices = [0 0 0; 0 4 0; 3 2 0;-3 2 0];
faces = [1 3 2;1 2 4];
mesh = GenericTriMesh(vertices, faces);
bnd = boundary(mesh);

assertEqual(testCase, 4, vertexNumber(bnd));
assertEqual(testCase, 4, edgeNumber(bnd));


function test_SimplyConnected(testCase) %#ok<*DEFNU>
% boundary of a simply connected mesh should have as many edges as vertices
ico = GenericTriMesh.createIcosahedron();
ico3s = smooth(subdivide(ico, 3));
box = [-0.5 2 -0.5 2 0.5 2];
icoC = clipVertices(ico3s, box);
bnd = boundary(icoC);

assertEqual(testCase, edgeNumber(bnd), vertexNumber(bnd));

