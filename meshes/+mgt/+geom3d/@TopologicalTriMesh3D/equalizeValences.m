function equalizeValences(obj)
%EQUALIZEVALENCES Equalize vertex valences by flipping edges.
%
%   output = equalizeValences(input)
%
%   Example
%   equalizeValences
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2023-04-06,    using Matlab 9.13.0.2049777 (R2022b)
% Copyright 2023 INRAE.

inds = find(obj.ValidEdges);

for i = 1:length(inds)
    ind = inds(i);

    faceInds = obj.EdgeFaces{ind};
    vertexInds = unique(obj.Faces(faceInds, :));
    
    % determine target valences of vertices (depends on whether vertices
    % are boundary or not)
    targetValences = [0 0 0 0];
    for iv = 1:4
        if isBoundaryVertex(obj, vertexInds(iv))
            targetValences(iv) = 4;
        else
            targetValences(iv) = 6;
        end
    end

    % compute valence deviation before flip
    valence0 = 0;
    for iv = 1:4
        deg = vertexDegree(obj, vertexInds(iv));
        valence0 = valence0 + abs(deg - targetValences(iv));
    end

    flipEdge(obj, ind);

    % compute valence deviation after flip
    valence1 = 0;
    for iv = 1:4
        deg = vertexDegree(obj, vertexInds(iv));
        valence1 = valence1 + abs(deg - targetValences(iv));
    end

    % if deviation was not improved, flip back
    if valence0 <= valence1
        flipEdge(obj, ind);
    end
end

