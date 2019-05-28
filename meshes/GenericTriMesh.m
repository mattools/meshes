classdef GenericTriMesh < handle
%GENERICTRIMESH A 3D triangular mesh with topology that can be modified 
%
%   Class GenericTriMesh
%
%   Example
%     % Create and display an octahedron
%     oct = GenericTriMesh.createOctahedron;
%     figure; hold on; axis equal; view(3)
%     draw(oct);
%
%     % Create and display an octahedron
%     ico = GenericTriMesh.createIcosahedron;
%     figure; hold on; axis equal; view(3)
%     draw(ico);
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-05-28,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The array of vertices, as a Nv-by-3 array of corodinates
    Vertices;
    
    % The topology on edges, as a NE-by-2 array containing vertex indices
    Edges;
    
    % the NF-by-3array containing vertex indices of each face.
    Faces;
    
    % cell array containing indices of faces associated to each edge.
    % should be 2 for regular edges, or 1 for boundary edges, but can be
    % more for non manifold edges.
    EdgeFaces = cell(1, 0);
    
    FaceEdges = cell(1, 0);
    
    % The list of edges adjacent to each vertex, as a 1-by-NV cell array
    VertexEdges = cell(1, 0);
    
    % The list of faces adjacent to each vertex, as a 1-by-NV cell array
    VertexFaces = cell(1, 0);
    
end % end properties

%% Static factories
methods (Static)
    function obj = createOctahedron()
        vertices = [1 0 0;0 1 0;-1 0 0;0 -1 0;0 0 1;0 0 -1];
        faces = [1 2 5;2 3 5;3 4 5;4 1 5;1 6 2;2 6 3;3 6 4;1 4 6];
        edges = [1 2;1 4;1 5; 1 6;2 3;2 5;2 6;3 4;3 5;3 6;4 5;4 6];
        obj = GenericTriMesh(vertices, faces);
        obj.Edges = edges;
    end
    
    function obj = createIcosahedron()
        len = 1/sin(pi/5)/2;
        z1 = sqrt(1-len*len);
        t1 = (0:2*pi/5:2*pi*(1-1/5))';
        x1 = len*cos(t1);  y1 = len*sin(t1);
        t2 = t1 + 2*pi/10;
        x2 = len*cos(t2); y2 = len*sin(t2);
        h = sqrt(len*len-.5*.5);
        z2 = sqrt(3/4 - (len-h)*(len-h));

        vertices = [0 0 0;...
            [x1 y1 repmat(z1, [5 1])]; ...
            [x2 y2 repmat(z1+z2, [5 1])]; ...
            0 0 2*z1+z2];
        edges = [...
            1 2;1 3;1 4;1 5;1 6; ...
            2 3;3 4;4 5;5 6;6 2; ...
            2 7;7 3;3 8;8 4;4 9;9 5;5 10;10 6;6 11;11 2; ...
            7 8;8 9;9 10;10 11;11 7; ...
            7 12;8 12;9 12;10 12;11 12];
        % faces are ordered to have normals pointing outside of the mesh
        faces = [...
            1 3  2 ; 1 4  3 ; 1  5  4 ;  1  6  5 ;  1 2  6;...
            2 3  7 ; 3 4  8 ; 4  5  9 ;  5  6 10 ;  6 2 11;...
            7 3  8 ; 8 4  9 ; 9  5 10 ; 10  6 11 ; 11 2  7;...
            7 8 12 ; 8 9 12 ; 9 10 12 ; 10 11 12 ; 11 7 12];

        obj = GenericTriMesh(vertices, faces);
        obj.Edges = edges;
    end
end % end methods


%% Constructor
methods
    function obj = GenericTriMesh(varargin)
    % Constructor for GenericTriMesh class

        if nargin == 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
            % classical vertex-face constructor
            obj.Vertices = varargin{1};
            obj.Faces = varargin{2};
            
        end
    end

end % end constructors

%% topology management
methods
    function [b1, b2] = isManifold(obj)
        %ISMANIFOLDMESH Check whether the mesh may be considered as manifold.
        %
        %   B = isManifoldMesh(MESH)
        %   Checks if the specified mesh is a manifold. When mesh is a manifold,
        %   all edges are connected to either 2 or 1 faces.
        %
        %   [B, HASBORDER] = isManifoldMesh(MESH)
        %   Also checks whether the mesh contains border faces. Border faces
        %   contains at least one edge which is ajacent to only one face.
        %
        %   Example
        %     mesh = GenericTriMesh.createOctahedron;
        %     isManifold(mesh)
        %     ans =
        %       logical
        %        1

        % make sure all inner data are up-to-date
        ensureValidEdges(obj);
        ensureValidFaceEdges(obj);
        ensureValidEdgeFaces(obj);
        
        % compute degree of each edge
        nEdges = edgeNumber(obj);
        edgeFaceDegrees = zeros(nEdges, 1);
        for iEdge = 1:nEdges
            edgeFaceDegrees(iEdge) = length(obj.EdgeFaces{iEdge});
        end
        
        % for each face, concatenate the face degree of each edge
        nFaces = faceNumber(obj);
        faceEdgeDegrees = zeros(nFaces, 3);
        for iFace = 1:nFaces
            edgeInds = obj.FaceEdges(iFace);
            faceEdgeDegrees(iFace, :) = edgeFaceDegrees(edgeInds);
        end
        
        regFaces = sum(ismember(faceEdgeDegrees, [1 2]), 2) == 3;
        innerFaces = sum(faceEdgeDegrees == 2, 2) == 3;
        borderFaces = regFaces & ~innerFaces;
        
        % check if mesh is manifold: all faces are either regular or border
        b1 = all(regFaces);
        
        % check if some faces are border
        b2 = any(borderFaces);
    end
end

%% Vertex management methods
methods
    function nv = vertexNumber(obj)
        nv = size(obj.Vertices, 1);
    end
end

%% Edge management methods
methods
    function ne = edgeNumber(obj)
        % ne = edgeNumber(mesh)
        if isempty(obj.Edges)
            computeEdgeVertices(obj);
        end
        ne = size(obj.Edges, 1);
    end
end

%% Face management methods
methods
    function nf = faceNumber(obj)
        nf = size(obj.Faces, 1);
    end
    
    function removeFace(obj, faceIndex)
        % Removes a face. This resets topological properties.
        obj.Faces(faceIndex, :) = [];
        clearEdges(obj);
    end
end


%% Local methods for updating local Properties 
methods
    function clearEdges(obj)
        % reset all data related to edges as well as VertexFaces.
        % (to be called when vertices of faces properties are updated)
        obj.Edges = zeros(0, 2);
        obj.EdgeFaces = cell(1, 0);
        obj.FaceEdges = zeros(0, 3);
        obj.VertexEdges = cell(1, 0);
        obj.VertexFaces = cell(1, 0);
    end
    
    function ensureValidEdges(obj)
        if isempty(obj.Edges)
            computeEdgeVertices(obj);
        end
    end
    
    function computeEdgeVertices(obj)
        % updates the property "Edges"
        
        % compute total number of edges
        % (3 edges per face)
        nFaces  = size(obj.Faces, 1);
        nEdges  = nFaces * 3;
        
        % create vertex indices for all edges (including duplicates)
        edges = zeros(nEdges, 2);
        for i = 1:nFaces
            f = obj.Faces(i, :);
            edges(((i-1)*3+1):i*3, :) = [f' f([2:end 1])'];
        end
        
        % remove duplicate edges, and sort the result
        obj.Edges = sortrows(unique(sort(edges, 2), 'rows'));
    end
    
    function ensureValidFaceEdges(obj)
        if isempty(obj.FaceEdges)
            computeFaceEdges(obj);
        end
    end
    
    function computeFaceEdges(obj)
        % compute face to edge indices array
        % as a Nf-by-3 array (each face connected to exactly three edges)
        
        nFaces = faceNumber(obj);
        obj.FaceEdges = zeros(nFaces, 3);
        
        % assume edges are sorted along dimension 2
        
        % iterate over faces to populate FaceEdges array
        for iFace = 1:nFaces
            % extract vertex indices of current face
            face = obj.Faces(iFace, :);
            
            % for each couple of adjacent vertices, find the index of the matching
            % row in the edges array
            for iEdge = 1:3
                % compute index of each edge vertex
                edge = sort([face(iEdge) face(mod(iEdge, 3) + 1)]);
                v1 = edge(1);
                v2 = edge(2);
                
                % find the matching row
                ind = find(obj.Edges(:,1) == v1 & obj.Edges(:,2) == v2);
                obj.FaceEdges(iFace, iEdge) = ind;
            end
        end
    end
    
    function ensureValidEdgeFaces(obj)
        if isempty(obj.EdgeFaces)
            computeEdgeFaces(obj);
        end
    end
    
    function computeEdgeFaces(obj)
        %Compute index of faces adjacent to each edge of a mesh.
        %
        %   Compute index array of faces adjacent to each edge of a mesh.
        %   The result EF has as many rows as the number of edges, and two column.
        %   The first column contains index of faces located on the left of the
        %   corresponding edge, whereas the second column contains index of the
        %   face located on the right. Some indices may be 0 if the mesh is not
        %   'closed'.
        
        % indices of faces adjacent to each edge
        obj.EdgeFaces = cell(1, edgeNumber(obj));
        
        for iFace = 1:faceNumber(obj)
            face = obj.Faces(iFace, :);
            
            % iterate on face edges
            for j = 1:length(face)
                % build edge: array of vertices
                j2 = mod(j, length(face)) + 1;
                
                % do not process edges with same vertices
                if face(j) == face(j2)
                    continue;
                end
                
                % vertex indices of current edge
                currentEdge = sort([face(j) face(j2)]);
                edgeIndex = find(ismember(obj.Edges, currentEdge, 'rows'));

                % check edge index could be identified
                if isempty(edgeIndex)
                    warning('GenericTriMesh:IllegalTopology', ...
                        'For face %d, edge %d (%d,%d) could not be found in edge array', ...
                        j, iFace, currentEdge(1), currentEdge());
                    continue;
                end
                
                % update list of face indices for current edge
                inds = obj.EdgeFaces{edgeIndex};
                if isempty(inds)
                    inds = iFace;
                else
                    inds = [inds iFace]; %#ok<AGROW>
                end
                obj.EdgeFaces{edgeIndex} = inds;
            end
        end
        
    end
end % end methods


%% Methods mimicking Geometry interface
methods
    function box = boundingBox(obj)
        % Returns the bounding box of this shape
        mini = min(obj.Vertices);
        maxi = max(obj.Vertices);
        box = Box3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function h = draw(varargin)
        % Draw the current geometry, eventually specifying the style
        
        % extract handle of axis to draw in
        if numel(varargin{1}) == 1 && ishghandle(varargin{1}, 'axes')
            hAx = varargin{1};
            varargin(1) = [];
        else
            hAx = gca;
        end

        % extract the point instance from the list of input arguments
        obj = varargin{1};
        varargin(1) = [];
        
        % add default drawing options
        options = {'FaceColor', [.75 .75 .75]};

        % extract optional drawing options
        if nargin > 1 && ischar(varargin{1})
            options = [options varargin];
        end
        
        if length(options) == 1
            options = [{'facecolor', [.75 .75 .75]} options];
        end

        h = patch('Parent', hAx, ...
            'vertices', obj.Vertices, 'faces', obj.Faces, ...
            options{:} );

        % optionnally add style processing
        if ~isempty(varargin) && isa(varargin{1}, 'Style')
            apply(varargin{1}, hh);
        end
                
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, varargin)
        % Returns a scaled version of this geometry
        factor = varargin{1};
        res = GenericTriMesh(obj.Vertices * factor, obj.Faces);
    end
    
    function res = translate(obj, varargin)
        % Returns a translated version of this geometry
        shift = varargin{1};
        res = GenericTriMesh(bsxfun(@plus, obj.Vertices, shift), obj.Faces);
    end
    
end % end methods

end % end classdef

