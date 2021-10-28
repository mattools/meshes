function res = subdivide(obj, n)
% Apply smoothing on mesh.
%
%   output = truc(input)
%
%   Example
%   truc
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-27,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

% compute the edge array
ensureValidEdges(obj);
ensureValidFaceEdges(obj);

% initialise the array of new vertices
vertices2 = obj.Vertices;

% several interpolated positions
t = linspace(0, 1, n + 1)';
coef2 = t(2:end-1);
coef1 = 1 - t(2:end-1);

% keep an array containing index of new vertices for each original edge
nEdges = edgeNumber(obj);
edgeNewVertexIndices = zeros(nEdges, n-1);

% create new vertices on each edge
for iEdge = 1:nEdges
    % extract each extremity as a point
    v1 = obj.Vertices(obj.Edges(iEdge, 1), :);
    v2 = obj.Vertices(obj.Edges(iEdge, 2), :);
    % compute new points
    newPoints = coef1 * v1 + coef2 * v2;
    % add new vertices, and keep their indices
    edgeNewVertexIndices(iEdge,:) = size(vertices2, 1) + (1:n-1);
    vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
end

faces2 = zeros(0, 3);
nFaces = faceNumber(obj);
for iFace = 1:nFaces
    % compute index of each corner vertex
    face = obj.Faces(iFace, :);
    iv1 = face(1);
    iv2 = face(2);
    iv3 = face(3);
    
    % compute index of each edge
    ie1 = obj.FaceEdges(iFace, 1);
    ie2 = obj.FaceEdges(iFace, 2);
    ie3 = obj.FaceEdges(iFace, 3);
    
    % indices of new vertices on edges
    edge1NewVertexIndices = edgeNewVertexIndices(ie1, :);
    edge2NewVertexIndices = edgeNewVertexIndices(ie2, :);
    edge3NewVertexIndices = edgeNewVertexIndices(ie3, :);
    
    % keep vertex 1 as reference for edges 1 and 3
    if obj.Edges(ie1, 1) ~= iv1
        edge1NewVertexIndices = edge1NewVertexIndices(end:-1:1);
    end
    if obj.Edges(ie3, 1) ~= iv1
        edge3NewVertexIndices = edge3NewVertexIndices(end:-1:1);
    end
    
    % create the first new face, on 'top' of the original face
    topVertexInds = [edge1NewVertexIndices(1) edge3NewVertexIndices(1)];
    newFace = [iv1 topVertexInds];
    faces2 = [faces2; newFace]; %#ok<AGROW>
    
    % iterate over middle strips
    for iStrip = 2:n-1
        % index of extreme vertices of current row
        ivr1 = edge1NewVertexIndices(iStrip);
        ivr2 = edge3NewVertexIndices(iStrip);
        
        % extreme vertices as points
        v1 = vertices2(ivr1, :);
        v2 = vertices2(ivr2, :);
        
        % create additional vertices within the bottom row of the strip
        t = linspace(0, 1, iStrip+1)';
        coef2 = t(2:end-1);
        coef1 = 1 - t(2:end-1);
        newPoints = coef1 * v1 + coef2 * v2;
        
        % compute indices of new vertices in result array
        newInds = size(vertices2, 1) + (1:iStrip-1);
        botVertexInds = [ivr1 newInds ivr2];
        
        % add new vertices
        vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
        
        % create top faces of current strip
        for k = 1:iStrip-1
            newFace = [topVertexInds(k) botVertexInds(k+1) topVertexInds(k+1)];
            faces2 = [faces2; newFace]; %#ok<AGROW>
        end
        
        % create bottom faces of current strip
        for k = 1:iStrip
            newFace = [topVertexInds(k) botVertexInds(k) botVertexInds(k+1)];
            faces2 = [faces2; newFace]; %#ok<AGROW>
        end
        
        % bottom vertices of current strip are top vertices of next strip
        topVertexInds = botVertexInds;
    end
    
    % for edge 2, keep vertex 2 of the current face as reference
    if obj.Edges(ie2, 1) ~= iv2
        edge2NewVertexIndices = edge2NewVertexIndices(end:-1:1);
    end
    
    % consider new vertices together with extremities
    botVertexInds = [iv2 edge2NewVertexIndices iv3];
    
    % create top faces for last strip
    for k = 1:n-1
        newFace = [topVertexInds(k) botVertexInds(k+1) topVertexInds(k+1)];
        faces2 = [faces2; newFace]; %#ok<AGROW>
    end
    
    % create bottom faces for last strip
    for k = 1:n
        newFace = [topVertexInds(k) botVertexInds(k) botVertexInds(k+1)];
        faces2 = [faces2; newFace]; %#ok<AGROW>
    end
end

res = TopologicalTriMesh3D(vertices2, faces2);