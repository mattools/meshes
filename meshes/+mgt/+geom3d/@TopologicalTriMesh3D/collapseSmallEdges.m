function collapseSmallEdges(obj, varargin)
%COLLAPSESMALLEDGES Collapse all edges smaller than a given value.
%
%   Usage:
%   collapseSmallEdges(MESH, MINLENGTH);
%   Update the mesh data structure by collapsing edges whose length is
%   smaller than value MINLENGTH.
%
%   collapseSmallEdges(MESH);
%   Automatically determines the minimum length of edges as a fraction of
%   the average length of all edges.
%
%   Example
%   collapseSmallEdges
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2023-04-06,    using Matlab 9.13.0.2049777 (R2022b)
% Copyright 2023 INRAE.

% determines the criterium used for selecting edges
[edgeLengths, inds] = edgeLength(obj);
if isempty(varargin)
    % use a fraction of the average length as length bound
    minLength = 0.10 * mean(edgeLengths);
elseif length(varargin) == 1
    % use numeric argument as length bound
    minLength = varargin{1};
else
    error('Unable to process input arguments');
end

% pre-select the edges to collapse
% (some of them may become invalid when a neighbor edge has collapsed)
indsToRemove = inds(edgeLengths < minLength);

% Iterate over edges
for iInd = 1:length(indsToRemove)
    edgeIndex = indsToRemove(iInd);
%     fprintf('iter %d, process edge#%d\n', iInd, edgeIndex);

    % process only valid edges
    % (selected edges may become invalid during merge process)
    if ~obj.ValidEdges(edgeIndex)
        continue;
    end

    % index of source and target vertices
    iv1 = obj.Edges(edgeIndex, 1);
    iv2 = obj.Edges(edgeIndex, 2);

    % case of boundary vertices
    if isBoundaryVertex(obj, iv1) && isBoundaryVertex(obj, iv2)
        if ~isBoundaryEdge(obj, edgeIndex)
            continue;
        end
    end

    % check that vertices adjacent to both extremities form an existing face
    % we can check only the faces incident to the edge to collapse
    neighs1 = unique(obj.Edges(obj.VertexEdges{iv1},:));
    neighs1(ismember(neighs1, [iv1 iv2])) = [];
    neighs2 = unique(obj.Edges(obj.VertexEdges{iv2},:));
    neighs2(ismember(neighs2, [iv1 iv2])) = [];
    adjBoth = neighs1(ismember(neighs1, neighs2));
    if length(adjBoth) > 2
        continue;
    end

    % call the collapse method
    collapseEdge(obj, edgeIndex);
end
