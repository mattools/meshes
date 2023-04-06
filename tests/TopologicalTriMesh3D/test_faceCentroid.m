function tests = test_faceCentroid
% Test suite for the file faceCentroid.
%
%   Test suite for the file faceCentroid
%
%   Example
%   test_faceCentroid
%
%   See also
%     faceCentroid

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>

% Test call of function without argument.
mesh = mgt.geom3d.TopologicalTriMesh3D();
addVertex(mesh, [0 0 0;12 0 0; 0 12 0]);
addFace(mesh, [1 2 3]);

centroids = faceCentroid(mesh);

assertEqual(testCase, size(centroids), [1 3]);
assertEqual(testCase, centroids, [4 4 0]);


function test_AfterRemovalOfFaces(testCase) %#ok<*DEFNU>

% Test call of function without argument.
mesh = mgt.geom3d.TopologicalTriMesh3D();
addVertex(mesh, [0 0 0;12 0 0; 0 12 0;12 12 0]);
addFace(mesh, [1 2 3]);
addFace(mesh, [3 2 4]);
removeFace(mesh, 1);

centroids = faceCentroid(mesh);

assertEqual(testCase, size(centroids), [1 3]);
assertEqual(testCase, centroids, [8 8 0]);

