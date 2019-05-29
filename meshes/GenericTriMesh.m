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
        % faces are ordered to have normals pointing outside of the mesh
        faces = [...
            1 3  2 ; 1 4  3 ; 1  5  4 ;  1  6  5 ;  1 2  6;...
            2 3  7 ; 3 4  8 ; 4  5  9 ;  5  6 10 ;  6 2 11;...
            7 3  8 ; 8 4  9 ; 9  5 10 ; 10  6 11 ; 11 2  7;...
            7 8 12 ; 8 9 12 ; 9 10 12 ; 10 11 12 ; 11 7 12];
        obj = GenericTriMesh(vertices, faces);
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


%% High-level processing
methods
   
    function res = smooth(obj, varargin)

        % determine number of iterations
        nIter = 1;
        if ~isempty(varargin)
            nIter = varargin{1};
        end

        % compute adjacency matrix,
        % result is a Nv-by-Nv matrix with zeros on the diagonal
        adj = adjacencyMatrix(obj);

        % ensure the size of the matrix is Nv-by-Nv
        % (this can not be the case if some vertices are not referenced)
        nv = vertexNumber(obj);
        if size(adj, 1) < nv
            adj(nv, nv) = 0;
        end

        % Add "self adjacencies"
        adj = adj + speye(nv);

        % weight each vertex by the number of its neighbors
        w = spdiags(full(sum(adj, 2).^(-1)), 0, nv, nv);
        adj = w * adj;

        % do averaging to smooth the field
        v2 = obj.Vertices;
        for k = 1:nIter
            v2 = adj * v2;
        end
        
        % create a new mesh with same faces
        res = GenericTriMesh(v2, obj.Faces);
    end

    function res = subdivide(obj, n)
        
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
        
        res = GenericTriMesh(vertices2, faces2);
    end
end


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
    
    function ind = addVertex(obj, position)
        if any(size(position) ~= [1 3])
            error('Require a 1-by-3 array of coordinates as input argument');
        end
        obj.Vertices = [obj.Vertices ; position];
        ind = size(obj.Vertices, 1);

        % optionnally updates topological data structures
        if ~isempty(obj.VertexEdges)
            obj.VertexEdges{ind} = [];
        end
        if ~isempty(obj.VertexFaces)
            obj.VertexFaces{ind} = [];
        end
    end
end

%% Edge management methods
methods
    function ne = edgeNumber(obj)
        % ne = edgeNumber(mesh)
        if isempty(obj.Edges)
            computeEdges(obj);
        end
        ne = size(obj.Edges, 1);
    end
    
    function ind = addEdge(obj, edgeVertices)
        % Adds an edge given vertex indices and return new edge index
        edge = sort(edgeVertices, 2); v1 = edge(1); v2 = edge(2);
        if ~isempty(obj.Edges)
            if ~isempty(find(obj.Edges(:,1) == v1 & obj.Edges(:,2) == v2, 1))
                error('Edge is already into mesh');
            end
        end
        obj.Edges = [obj.Edges; edge];
        ind = size(obj.Edges, 1);
        
        % optionnally add edge index to vertices
        if ~isempty(obj.VertexEdges)
            obj.VertexEdges{v1} = [obj.VertexEdges{v1} ind];
            obj.VertexEdges{v2} = [obj.VertexEdges{v2} ind];
        end
        
        % optionnally create empty list of faces for this edge
        if ~isempty(obj.EdgeFaces)
            obj.EdgeFaces{ind} = [];
        end
    end
end

%% Face management methods
methods
    function nf = faceNumber(obj)
        nf = size(obj.Faces, 1);
    end
    
    function indFace = addFace(obj, vertexInds)
        % Adds a face given 3 vertex indices and return new face index
        
        if any(size(vertexInds) ~= [1 3])
            error('Require an array of three indices as input argument');
        end
        obj.Faces = [obj.Faces ; vertexInds];
        indFace = size(obj.Faces, 1);

        % optionnaly updates vertex faces
        if ~isempty(obj.VertexFaces)
            obj.VertexFaces{face(1)} = [obj.VertexFaces{face(1)} ind];
            obj.VertexFaces{face(2)} = [obj.VertexFaces{face(2)} ind];
            obj.VertexFaces{face(3)} = [obj.VertexFaces{face(3)} ind];
        end
        
        % optionnaly updates edges
        if ~isempty(obj.Edges)
            newEdges = sort(vertexInds([1 2; 2 3;1 3]), 2);
            edgeInds = zeros(1,3);
            
            for iEdge = 1:3
                v1 = newEdges(iEdge, 1);
                v2 = newEdges(iEdge, 2);
                edgeInd = find(obj.Edges(:,1) == v1 & obj.Edges(:,2) == v2);
                if isempty(edgeInd)
                    edgeInd = addEdge(obj, [v1 v2]);
                end
                edgeInds(iEdge) = edgeInd;
            end
            
            if ~isempty(obj.FaceEdges)
                obj.FaceEdges(indFace, :) = edgeInds;
            end
            if ~isempty(obj.EdgeFaces)
                obj.EdgeFaces{edgeInds(1)} = [obj.EdgeFaces{edgeInds(1)} ind];
                obj.EdgeFaces{edgeInds(2)} = [obj.EdgeFaces{edgeInds(2)} ind];
                obj.EdgeFaces{edgeInds(3)} = [obj.EdgeFaces{edgeInds(3)} ind];
            end
        end
    end
    
    function removeFace(obj, faceIndex)
        % Removes a face. This resets topological properties.
        obj.Faces(faceIndex, :) = [];
        
        % TODO: update edges instead of clearing
        clearEdges(obj);
    end
    
    function adj = adjacencyMatrix(obj)

        % forces faces to be floating point array, for sparse function
        faces = obj.Faces;
        if ~isfloat(faces)
            faces = double(faces);
        end
        
        % populate a sparse matrix
        adj = sparse(...
            [faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
            [faces(:,3); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1)], ...
            1.0);
        
        % remove double adjacencies
        adj = min(adj, 1);
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
            computeEdges(obj);
        end
    end
    
    function computeEdges(obj)
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

