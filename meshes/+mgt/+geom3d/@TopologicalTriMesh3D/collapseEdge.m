function collapseEdge(obj, edgeIndex)
% Collapse the two extremities of the edge, removing adjacent faces.
%
%   collapseEdge(MESH, EDGE_IDX);
%   Collapses the two extremities of the edge specified by index EDGE_IDX,
%   and replaces the new vertex at the centroid of the two previous
%   vertices.
%   This removes one vertex, three edges, and two faces from the original
%   mesh.
%
%   Two conditions have to be met:
%   * If the extremity vertices are boundary vertices, then the edge must
%       be a boundary edge
%   * for all vertices incident to both extremity vertices, there must be a
%       face with the two extremity vertices and the incident vertex. If
%       not verified, this creates a dupicate fold-over triangle.
%
%   Example
%   
%
%   See also
%     removeFace

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.


%% Retrieve necessary info

% list of edges to remove, initializes with edge index.
edgesToRemove = edgeIndex;

% index of source and target vertices
iv1 = obj.Edges(edgeIndex, 1);
iv2 = obj.Edges(edgeIndex, 2);

% get index of faces adjacent to current edge
adjFaceInds = obj.EdgeFaces{edgeIndex};
if length(adjFaceInds) > 2
    warning('edge %d was adjacent to more than two faces', edgeIndex);
end


%% check pre-conditions

% case of boundary vertices
if isBoundaryVertex(obj, iv1) && isBoundaryVertex(obj, iv2)
    if ~isBoundaryEdge(obj, edgeIndex)
        error('Edge to remove (#%d) must be a boundary edge', edgeIndex);
    end
end

% check vertices adjacent to both extremities form an existing face
% we can check only the faces incident to the edge to collapse
edgeFaceVerts = obj.Faces(obj.EdgeFaces{edgeIndex},:);
neighs1 = unique(obj.Edges(obj.VertexEdges{iv1},:));
neighs1(ismember(neighs1, [iv1 iv2])) = [];
neighs2 = unique(obj.Edges(obj.VertexEdges{iv2},:));
neighs2(ismember(neighs2, [iv1 iv2])) = [];
adjBoth = neighs1(ismember(neighs1, neighs2));
for iAdj = 1:length(adjBoth)
    ivAdj = adjBoth(iAdj);
    face = [iv1 iv2 ivAdj];
    if ~any(sum(ismember(edgeFaceVerts, face), 2) == 3)
        error('Vertex #%d is adjacent to both extremities of edge #%d, but no face exist with all three vertices.', ivAdj, edgeIndex);
    end
end


%% Iterate over adjacent faces

% for each adjacent face:
% * remove ref from adjacent vertices
% * merge remaining edges
% * for face incident to edge2, update edge refs
% * invalidate face
for iAdjFace = 1:length(adjFaceInds)
    iFace = adjFaceInds(iAdjFace);
    % index of the vertex not belonging to the edge to collapse
    faceVertexInds = obj.Faces(iFace, :);
    iv3 = faceVertexInds(~ismember(faceVertexInds, [iv1 iv2]));
    
    % for each vertex incident to current face, remove ref to current face
    obj.VertexFaces{iv1}(obj.VertexFaces{iv1} == iFace) = [];
    obj.VertexFaces{iv3}(obj.VertexFaces{iv3} == iFace) = [];
    
    % index of the edges to keep and to remove
    ind1 = findEdgeIndex(obj, [iv1 iv3]); % to keep
    ind2 = findEdgeIndex(obj, [iv2 iv3]); % to remove
    if length(ind2) ~= 1
        error('Should have only one edge to remove!');
    end
    edgesToRemove = [edgesToRemove ; ind2]; %#ok<AGROW>
    
    % remove ref to edge2 from vertex v3
    obj.VertexEdges{iv3}(obj.VertexEdges{iv3} == ind2) = [];

    % for each face incident to edge ind2, need to update ref to edge ind1
    faceInds = obj.EdgeFaces{ind2};
    faceInds(faceInds == iFace) = [];
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
    faceInds(faceInds == iFace) = [];
    obj.EdgeFaces{ind1} = faceInds;
    obj.EdgeFaces{ind2} = [];
    
    % invalidate current face
    obj.FaceEdges(iFace, :) = [0 0 0];
    obj.ValidFaces(iFace) = false;
end


%% Clean-up 

% as v2 becomes v1, need to update refs to v2.
v2EdgeArray = obj.Edges(obj.VertexEdges{iv2}, :);
v2EdgeArray(v2EdgeArray == iv2) = iv1;
obj.Edges(obj.VertexEdges{iv2}, :) = v2EdgeArray;
v2FaceArray = obj.Faces(obj.VertexFaces{iv2}, :);
v2FaceArray(v2FaceArray == iv2) = iv1;
obj.Faces(obj.VertexFaces{iv2}, :) = v2FaceArray;

% update mapping from vertices to incident faces
inds = [obj.VertexFaces{iv1} obj.VertexFaces{iv2}];
inds(ismember(inds, adjFaceInds)) = [];
obj.VertexFaces{iv1} = inds;
obj.VertexFaces{iv2} = [];

% update mapping from vertices to incident edges
inds = [obj.VertexEdges{iv1} obj.VertexEdges{iv2}];
inds(ismember(inds, edgesToRemove)) = [];
obj.VertexEdges{iv1} = inds;
obj.VertexEdges{iv2} = [];

% remove old edges
obj.EdgeFaces{edgeIndex} = [];
obj.ValidEdges(edgesToRemove) = false;

% remove vertex v2
obj.ValidVertices(iv2) = false;
