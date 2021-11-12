function res = planeIntersection(obj, plane, varargin)
% Intersection of this mesh with a plane.
%
%   output = planeIntersection(MESH, PLANE)
%   Computes the intersection between a plane and a mesh given by vertex and
%   face lists. The result is an array of LinearRing3D objects.
%
%   Example
%     % Intersect a simple polyhedron with a plane
%     mesh = TopologicalTriMesh3D.createIcosahedron;
%     plane = Plane3D(Point3D([0.5 0.5 0.5]), Vector3D([3 4 5]));
%     % draw the primitives
%     figure; hold on; set(gcf, 'renderer', 'opengl');
%     axis equal; axis([-1 1 -1 1 0 2]); view(3);
%     draw(mesh); draw(plane);
%     % compute intersection polygon
%     intersection = planeIntersection(mesh, plane);
%     draw(intersection, 'LineWidth', 2);
%
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-11-10,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.


% identify which edges cross the mesh
% (one vertex is below, the other one is above)
inds = isBelow(plane, Point3D(obj.Vertices));
crossEdgeInds = find(sum(inds(obj.Edges), 2) == 1);
crossEdgeVerts = obj.Edges(crossEdgeInds, :);

% compute one intersection point for each edge
edges = LineSegment3D(obj.Vertices(crossEdgeVerts(:, 1), :), obj.Vertices(crossEdgeVerts(:, 2), :));
intersectionPoints = planeIntersection(edges, plane);

nCrossEdges = length(crossEdgeInds);

% create mapping from mesh edge indices to intersection indices
edgeIntersectionMap = containers.Map('keyType', 'int32', 'valueType', 'int32');
for iInter = 1:nCrossEdges
    edgeIntersectionMap(crossEdgeInds(iInter)) = iInter;
end

% create empty cell array of polygons
polys = {};


%% Iterate on edges and faces to form polygons

% initialize an array indicating which indices need to be processed
remainingCrossEdges = true(nCrossEdges, 1);

% iterate while there are some crossing edges to process
% while ~isempty(crossEdgesIndsToProcess)
while any(remainingCrossEdges)
    
    % start at any edge, mark it as current
    startCrossEdgeIndex = find(remainingCrossEdges, 1, 'first');
    % mark current edge as processed
    remainingCrossEdges(startCrossEdgeIndex) = false;
    
    % initialize new set of edge indices
    polyCrossEdgeInds = startCrossEdgeIndex;
    
    % convert to mesh edge index
    startEdgeIndex = crossEdgeInds(startCrossEdgeIndex);
    currentEdgeIndex = startEdgeIndex;
    
    % choose one of the two faces around the edge
    currentFaceIndex = obj.EdgeFaces{currentEdgeIndex}(1);

    % iterate along current face-edge couples until back to first edge
    while true
        % find the index of next crossing edge
        inds = obj.FaceEdges(currentFaceIndex, :);
        inds = inds(ismember(inds, crossEdgeInds));
        currentEdgeIndex = inds(inds ~= currentEdgeIndex);
     
        % mark current edge as processed
        currentCrossEdgeIndex = edgeIntersectionMap(currentEdgeIndex);
        remainingCrossEdges(currentCrossEdgeIndex) = false;
        
        % check end of current loop
        if currentEdgeIndex == startEdgeIndex
            break;
        end
        
        % find the index of the other face containing current edge
        inds = obj.EdgeFaces{currentEdgeIndex};
        currentFaceIndex = inds(inds ~= currentFaceIndex);
        
        % add index of current edge
        polyCrossEdgeInds = [polyCrossEdgeInds currentCrossEdgeIndex]; %#ok<AGROW>
    end
    
    % create polygon, and add it to list of polygons
    poly = intersectionPoints(polyCrossEdgeInds, :);
    polys = [polys, {poly}]; %#ok<AGROW>
end

% format output as an array of LinearRing3D objects
nPolys = length(polys);
res(nPolys, 1) = LinearRing3D();
for i = 1:length(polys)
    res(i) = LinearRing3D(polys{i});
end
