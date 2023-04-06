function ccList = connectedComponents(obj)
% Set of edge-connected components of the mesh.
%
%   CCList = connectedComponents(MESH)
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
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.


ensureValidFaceEdges(obj);
ensureValidEdgeFaces(obj);

nFaces = size(obj.Faces, 1);

% associate a component label to each face.
%  0 -> not yet assignd
% -1 -> invalid face (no vertex associated)
faceLabels = zeros(nFaces, 1);

% init
currentLabel = 0;
ccList = {};

while true
    % find a face without label
    facesToUpdate = find(faceLabels == 0, 1);
    
    if isempty(facesToUpdate)
        break;
    end
    
    currentLabel = currentLabel + 1;
    
    while ~isempty(facesToUpdate)
        indFace = facesToUpdate(1);
        facesToUpdate(1) = [];
        
        faceLabels(indFace) = currentLabel;
        
        edgeInds = obj.FaceEdges(indFace, :);
        for iEdge = 1:length(edgeInds)
            faceInds = obj.EdgeFaces{edgeInds(iEdge)};
            faceInds(faceInds == indFace) = [];
            
            % length of faceInds should be 0 or 1 for manifold meshes
            % but we keep the loop to manage non-manifold cases
            for iFace = 1:length(faceInds)
                if faceLabels(faceInds(iFace)) == 0
                    facesToUpdate = [facesToUpdate faceInds(iFace)]; %#ok<AGROW>
                end
            end
        end
    end
    
    % create new mesh with only necessary vertices and faces
    newFaces = obj.Faces(faceLabels == currentLabel, :);
    cc = trimmedMesh(TopologicalTriMesh3D(obj.Vertices, newFaces));
    ccList = [ccList {cc}]; %#ok<AGROW>
end
