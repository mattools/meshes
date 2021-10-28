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

% check input
if isempty(obj.Edges)
    error('Expect edge array to have been computed');
end
ensureValidEdgeFaces(obj);

% list of edges to remove
edgesToRemove = edgeIndex;

% index of source and target vertices
v1 = obj.Edges(edgeIndex, 1);
v2 = obj.Edges(edgeIndex, 2);

% get index of faces adjacent to current edge
adjFaceInds = obj.EdgeFaces{edgeIndex};

% for each adjacent face, merge remaining edges
for iFace = 1:length(adjFaceInds)
    currentFaceIdx = adjFaceInds(iFace);
    % index of the vertex not belonging to the edge to collapse
    faceVertexInds = obj.Faces(currentFaceIdx, :);
    v3 = faceVertexInds(~ismember(faceVertexInds, [v1 v2]));
    
    % index of the edges to keep and to remove
    ind1 = find(sum(ismember(obj.Edges, [v1 v3]), 2) == 2); % to keep
    ind2 = find(sum(ismember(obj.Edges, [v2 v3]), 2) == 2); % to remove
    
    if length(ind2) ~= 1
        error('Should have only one edge to remove!');
    end
    
    % update the faces incident to the edge to remove
    % (if the "FaceEdges" property was already defined)
    if ~isempty(obj.FaceEdges)
        % clear info for current face
        obj.FaceEdges(currentFaceIdx) = -1;
        
        % update the faces incident to the edge to remove
        faceInds = obj.EdgeFaces{ind2};
        faceInds(faceInds == currentFaceIdx) = [];
        for iFace2 = 1:length(faceInds) % should have only 1
            edgeInds2 = obj.FaceEdges(faceInds(iFace2), :);
            edgeInds2(edgeInds2 == ind2) = ind1;
            obj.FaceEdges(faceInds(iFace2), :) = edgeInds2;
        end
    end
    
    % update EdgeFaces info, by merging indices of all faces adjacent to
    % either edge1 or edge2 (excluding current face) 
    faceInds = [obj.EdgeFaces{ind1} obj.EdgeFaces{ind2}];
    faceInds(faceInds == currentFaceIdx) = [];
    obj.EdgeFaces{ind1} = faceInds;
    obj.EdgeFaces{ind2} = [];
end

% as v2 becomes v1, need to update refs to v2.
obj.Edges(obj.Edges == v2) = v1;
obj.Faces(obj.Faces == v2) = v1;

% remove vertex v2
obj.VertexValidities(v2) = false;

% update mapping from vertices to incident edges
if ~isempty(obj.VertexEdges)
    inds = [obj.VertexEdges{v1} obj.VertexEdges{v2}];
    inds(ismember(inds, edgesToRemove)) = [];
    obj.VertexEdges{v1} = inds;
    obj.VertexEdges{v2} = [];
end

% update mapping from vertices to incident faces
if ~isempty(obj.VertexFaces)
    inds = [obj.VertexFaces{v1} obj.VertexFaces{v2}];
    inds(ismember(inds, adjFaceInds)) = [];
    obj.VertexFaces{v1} = inds;
    obj.VertexFaces{v2} = [];
end

% mark vertices of removed faces as 'removed'
obj.Faces(adjFaceInds, :) = -1;

% remove old edges
obj.Edges(edgesToRemove, :) = -1;
