function collapseEdge(obj, edgeIndex)
% Collapse the two extremities of the edge, removing adjacent faces.
%
%   collapseEdge(MESH, EIDX);
%
%   Example
%   tre
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

% list of edges to remove
edgesToRemove = edgeIndex;

% index of source and target vertices
v1 = obj.Edges(edgeIndex, 1);
v2 = obj.Edges(edgeIndex, 2);

% get index of faces adjacent to current edge
adjFaceInds = obj.EdgeFaces{edgeIndex};
if length(adjFaceInds) > 2
    warning('edge %d was adjacent to more than two faces', edgeIndex);
end

% for each adjacent face:
% * remove ref from adjacent vertices
% * merge remaining edges
% * for face incident to edge2, update edge refs
% * invalidate face
for iFace = 1:length(adjFaceInds)
    currentFaceIdx = adjFaceInds(iFace);
    % index of the vertex not belonging to the edge to collapse
    faceVertexInds = obj.Faces(currentFaceIdx, :);
    v3 = faceVertexInds(~ismember(faceVertexInds, [v1 v2]));
    
    % for each vertex incident to current face, remove ref to current face
    obj.VertexFaces{v1}(obj.VertexFaces{v1} == currentFaceIdx) = [];
    obj.VertexFaces{v3}(obj.VertexFaces{v3} == currentFaceIdx) = [];
    
    % index of the edges to keep and to remove
    ind1 = findEdgeIndex(obj, [v1 v3]); % to keep
    ind2 = findEdgeIndex(obj, [v2 v3]); % to remove
    if length(ind2) ~= 1
        error('Should have only one edge to remove!');
    end
    edgesToRemove = [edgesToRemove ; ind2]; %#ok<AGROW>
    
    % remove ref to edge2 from vertex v3
    obj.VertexEdges{v3}(obj.VertexEdges{v3} == ind2) = [];

    % for each face incident to edge ind2, need to update ref to edge ind1
    faceInds = obj.EdgeFaces{ind2};
    faceInds(faceInds == currentFaceIdx) = [];
    if length(faceInds) > 1
        warning('edge %d was adjacent to more than two faces', ind2);
    end
    for iFace2 = 1:length(faceInds) % should have only 1
        edgeInds2 = obj.FaceEdges(faceInds(iFace2), :);
        edgeInds2(edgeInds2 == ind2) = ind1;
        obj.FaceEdges(faceInds(iFace2), :) = edgeInds2;
    end
    
    % update EdgeFaces info, by merging indices of all faces adjacent to
    % either edge1 or edge2 (excluding current face) 
    faceInds = [obj.EdgeFaces{ind1} obj.EdgeFaces{ind2}];
    faceInds(faceInds == currentFaceIdx) = [];
    obj.EdgeFaces{ind1} = faceInds;
    obj.EdgeFaces{ind2} = [];
    
    % invalidate current face
    obj.FaceEdges(currentFaceIdx, :) = [0 0 0];
    obj.ValidFaces(currentFaceIdx) = false;
end

% as v2 becomes v1, need to update refs to v2.
v2EdgeArray = obj.Edges(obj.VertexEdges{v2}, :);
v2EdgeArray(v2EdgeArray == v2) = v1;
obj.Edges(obj.VertexEdges{v2}, :) = v2EdgeArray;
v2FaceArray = obj.Faces(obj.VertexFaces{v2}, :);
v2FaceArray(v2FaceArray == v2) = v1;
obj.Faces(obj.VertexFaces{v2}, :) = v2FaceArray;

% update mapping from vertices to incident faces
inds = [obj.VertexFaces{v1} obj.VertexFaces{v2}];
inds(ismember(inds, adjFaceInds)) = [];
obj.VertexFaces{v1} = inds;
obj.VertexFaces{v2} = [];

% update mapping from vertices to incident edges
inds = [obj.VertexEdges{v1} obj.VertexEdges{v2}];
inds(ismember(inds, edgesToRemove)) = [];
obj.VertexEdges{v1} = inds;
obj.VertexEdges{v2} = [];

% remove old edges
obj.EdgeFaces{edgeIndex} = [];
obj.ValidEdges(edgesToRemove) = false;

% remove vertex v2
obj.ValidVertices(v2) = false;
