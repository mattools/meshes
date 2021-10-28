classdef TopologicalTriMesh3D < handle
% A 3D triangular mesh that stores the full topology information.
%
%   The TopologicalTriMesh3D class is an implementation of 3D triangular
%   mesh that stores the full topology information between vertices, edges
%   and faces. It should be possible to retrieve with edges are incident to
%   a vertex, which faces are incident to an edge, and so on...
%   It is expected to be possible to remove elements from the mesh (work in
%   progress).
%
%   Requires the "Geometry" package.
%
%   Example
%     % Create and display an octahedron
%     oct = TopologicalTriMesh3D.createOctahedron;
%     figure; hold on; axis equal; view(3)
%     draw(oct);
%
%     % Create and display an octahedron
%     ico = TopologicalTriMesh3D.createIcosahedron;
%     figure; hold on; axis equal; view(3)
%     draw(ico);
%
%   See also
%     TriMesh3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-05-28,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The array of vertices, as a Nv-by-3 array of coordinates.
    Vertices;
    
    % An array of boolean flags indicating for which vertex if it is valid.
    % (default is true for each vertex)
    VertexValidities;
    
    % The topology on edges, as a NE-by-2 array containing vertex indices.
    Edges;
    
    % The NF-by-3 array containing vertex indices of each face.
    % Some rows may contains [-1 -1 -1] if the face was removed.
    Faces;
    
    % Cell array containing indices of faces associated to each edge.
    % should be 2 for regular edges, or 1 for boundary edges, but can be
    % more for non manifold edges.
    EdgeFaces = cell(1, 0);
    
    FaceEdges = cell(1, 0);
    
    % The list of edges adjacent to each vertex, as a 1-by-NV cell array.
    VertexEdges = cell(1, 0);
    
    % The list of faces adjacent to each vertex, as a 1-by-NV cell array.
    VertexFaces = cell(1, 0);
    
end % end properties

%% Static factories
methods (Static)
    function obj = createOctahedron()
        vertices = [1 0 0;0 1 0;-1 0 0;0 -1 0;0 0 1;0 0 -1];
        faces = [1 2 5;2 3 5;3 4 5;4 1 5;1 6 2;2 6 3;3 6 4;1 4 6];
        edges = [1 2;1 4;1 5; 1 6;2 3;2 5;2 6;3 4;3 5;3 6;4 5;4 6];
        obj = TopologicalTriMesh3D(vertices, faces);
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
        obj = TopologicalTriMesh3D(vertices, faces);
    end
end % end methods


%% Constructor
methods
    function obj = TopologicalTriMesh3D(varargin)
        % Constructor for TopologicalTriMesh3D class.

        if nargin == 0
            % empty constructor -> nothing to do
            
        elseif nargin == 1 && isa(varargin{1}, 'TriMesh3D')
            % convert TriMesh to topological mesh
            var1 = varargin{1};
            setVertices(obj, var1.Vertices);
            obj.Faces = var1.Faces;
            
        elseif nargin == 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
            % classical vertex-face constructor
            setVertices(obj, varargin{1});
            obj.Faces = varargin{2};
        else
            error('Unable to interpret input arguments');
            
        end
    end

end % end constructors


%% High-level processing
methods
    function res = reverseOrientation(obj)
        % Reverse the orientation of the normals of the mesh.
        %
        %    MESH2 = reverseOrientation(MESH);
        faces2 = obj.Faces(:, [1 3 2]);
        res = TopologicalTriMesh3D(obj.Vertices, faces2);
    end
end


%% Mesh edition methods
methods
    function splitEdge(obj, edgeIndex)
        
        % get indices of faces adjacent to current edge
        ensureValidEdgeFaces(obj);
        adjFaceInds = obj.EdgeFaces{edgeIndex};
        
        % check local topology
        if length(adjFaceInds) ~= 2
            error('Can only split edges with two adjacent faces, not %d', length(adjFaceInds));
        end
        f1 = adjFaceInds(1);
        f2 = adjFaceInds(2);

        % get indices of adjacent and opoosite vertices around current edge
        v1 = obj.Edges(edgeIndex, 1);
        v2 = obj.Edges(edgeIndex, 2);
        face1 = obj.Faces(f1, :);
        vf1 = face1(~ismember(face1, [v1 v2]));
        face2 = obj.Faces(f2, :);
        vf2 = face2(~ismember(face2, [v1 v2]));
        
        % add new vertex in the middle of new edge
        pos = edgeCentroid(obj, edgeIndex);
        vNew = addVertex(obj, pos);
        
        % update first face: replace v2 with vNew
        face1(face1 == v2) = vNew;
        obj.Faces(f1, :) = face1;
        % update second face: replace v2 with vNew
        face2(face2 == v2) = vNew;
        obj.Faces(f2, :) = face2;
        
        % add three new edges
        ne0 = addEdge(obj, [vNew v2]);
        ne1 = addEdge(obj, [vNew vf1]);
        ne2 = addEdge(obj, [vNew vf2]);
        
        % add three new faces
        nf1 = addFace(obj, [vf1 v2 vNew]);
        nf2 = addFace(obj, [vf2 v2 vNew]);
        
        %              vf1                               vf1
        %           -       -                         -   |   -
        %       e11     f1    e21                 e11    ne1    e21
        %     -                   -             -    f1   |  nf1    -
        % v1 ----------------------- v2     v1 -- e0 -- vNew -- ne0 -- v2
        %     -                   -             -    f2   |  nf2    -
        %       e12     f2    e22                 e12    ne2    e22
        %           -       -                         -   |   -
        %              vf2                               vf2
        
        % attach new faces to new edges
        obj.EdgeFaces{ne0} = [nf1, nf2];
        obj.EdgeFaces{ne1} = [f1, nf1];
        obj.EdgeFaces{ne2} = [f2, nf2];
        
        % TODO: if necessary, updates faceEdges info
        % TODO: if necessary, updates vertexEdges info
    end
    
    function splitFace(obj, faceIndex)
        % Splits the specified face into three triangular faces.
        centro = faceCentroids(obj, faceIndex);
        vNew = addVertex(obj, centro);
        v1 = obj.Faces(faceIndex, 1);
        v2 = obj.Faces(faceIndex, 2);
        v3 = obj.Faces(faceIndex, 3);
        obj.Faces(faceIndex, 3) = vNew;
        addFace(obj, [v2 v3 vNew]);
        addFace(obj, [v3 v1 vNew]);
        
        clearEdges(obj);
        
        % TODO: if necessary, updates Edges info
        % TODO: if necessary, updates faceEdges info
        % TODO: if necessary, updates vertexEdges info
    end
    
    function flipEdge(obj, edgeIndex)
        %
        %              vf1                               vf1
        %           -       -                         -   |   -
        %       e11     f1    e21                 e11     |     e21
        %     -                   -             -         |         -
        % v1 ---------- e ---------- v2     v1      f1    e    f2      v2
        %     -                   -             -         |         -
        %       e12     f2    e22                 e12     |     e22
        %           -       -                         -   |   -
        %              vf2                               vf2
        
        % get indices of faces adjacent to current edge
        ensureValidEdgeFaces(obj);
        adjFaceInds = obj.EdgeFaces{edgeIndex};
        
        % check local topology
        if length(adjFaceInds) ~= 2
            error('Can only flip edges with two adjacent faces, not %d', length(adjFaceInds));
        end
        f1 = adjFaceInds(1);
        f2 = adjFaceInds(2);

        % get indices of adjacent and opoosite vertices around current edge
        v1 = obj.Edges(edgeIndex, 1);
        v2 = obj.Edges(edgeIndex, 2);
        face1 = obj.Faces(f1, :);
        vf1 = face1(~ismember(face1, [v1 v2]));
        face2 = obj.Faces(f2, :);
        vf2 = face2(~ismember(face2, [v1 v2]));

        obj.Faces(f1, :) = [v1 vf2 vf1];
        obj.Faces(f2, :) = [v2 vf1 vf2];

        obj.Edges(edgeIndex, :) = sort([vf1 vf2]);
        
        % TODO: if necessary, updates faceEdges info
        % TODO: if necessary, updates vertexEdges info
        
    end
end

%% Geometry methods
methods
    function box = bounds(obj)
        % Return the bounds of this mesh.
        mini = min(obj.Vertices);
        maxi = max(obj.Vertices);
        box = Bounds3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function lengths = edgeLength(obj, varargin)
        ensureValidEdges(obj);
        lengths = sum((obj.Vertices(obj.Edges(:,1), :) - obj.Vertices(obj.Edges(:,2), :)).^2, 2);
    end
    
    function centroids = faceCentroids(obj, varargin)
        % Compute the normals of all faces in the mesh
        
        if nargin == 1
            nf = size(obj.Faces, 1);
            centroids = zeros(nf, 3);
            % For triangular meshes, uses accelerated method
            % (taken from https://github.com/alecjacobson/gptoolbox)
            for ff = 1:3
                centroids = centroids + obj.Vertices(obj.Faces(:,ff),:) / 3.0;
            end
        else
            indFaces = varargin{1};
            nf = length(indFaces);
            centroids = zeros(nf, 3);
            % For triangular meshes, uses accelerated method
            % (taken from https://github.com/alecjacobson/gptoolbox)
            for ff = 1:3
                centroids = centroids + obj.Vertices(obj.Faces(indFaces,ff),:) / 3.0;
            end
        end
    end
    
    function normals = faceNormals(obj)
        % Compute the normals of faces.

        % compute vector of first edges
        v1 = obj.Vertices(obj.Faces(:,1),:);
        v12 = obj.Vertices(obj.Faces(:,2),:) - v1;
        v13 = obj.Vertices(obj.Faces(:,3),:) - v1;
        
        % compute normals using cross product (vectors have same size)
        normals = cross(v12, v13, 2);
    end
    
    function normals = vertexNormals(obj)
        % Compute the normals at vertices.
        nv = size(obj.Vertices, 1);
        nf = size(obj.Faces, 1);
        
        % unit normals to the faces
        fNormals = faceNormals(obj);
        
        % compute normal of each vertex: sum of normals to each face
        normals = zeros(nv, 3);
        for i = 1:nf
            face = obj.Faces(i, :);
            for j = 1:length(face)
                v = face(j);
                normals(v, :) = normals(v,:) + fNormals(i,:);
            end
        end
    end

    function centro = edgeCentroid(obj, edgeIndex)
        % Compute centroid of edge.
        
        ensureValidEdges(obj);
        v1 = obj.Vertices(obj.Edges(edgeIndex,1), :);
        v2 = obj.Vertices(obj.Edges(edgeIndex,2), :);
        centro = (v1 + v2) / 2;
    end
end


%% Topological queries
methods
    function res = vertexLink(obj, vertexIndex)
        % Return the 1-link around the specifid vertex as a new mesh.
        
        ensureValidEdges(obj);
        edgeInds = sum(obj.Edges == vertexIndex, 2) > 0;
        vertexInds = unique(obj.Edges(edgeInds, :));
        vertexInds(vertexInds == vertexIndex) = [];
        
        newVertices = obj.Vertices(vertexInds, :);
        
        % find index of edges belonging to the new complex
        linkEdgeInds = sum(ismember(obj.Edges, vertexInds), 2) == 2;
        newEdges = obj.Edges(linkEdgeInds, :);
        for i = 1:numel(newEdges)
            newEdges(i) = find(vertexInds == newEdges(i), 1);
        end
        
        % create the resulting mesh
        res = TopologicalTriMesh3D();
        setVertices(res, newVertices);
        res.Edges = newEdges;
    end
    
    function polyList = vertexLinkPolygons(obj, vertexIndex)
        % Return the 1-link around a vertex as a list of 3D polylines.
        linkMesh = vertexLink(obj, vertexIndex);

        edges = linkMesh.Edges;
        polyList = cell(0, 0);
        while ~isempty(edges)
            v0 = edges(1,1);
            v1 = edges(1,2);
            edges(1,:) = [];

            poly = linkMesh.Vertices(v0, :);
            while v1 ~= v0
                poly = [poly ; linkMesh.Vertices(v1, :)]; %#ok<AGROW>
                indEdge = sum(edges == v1, 2) > 0;
                edge = edges(indEdge, :);
                edges(indEdge,:) = [];
                v1 = edge(edge ~= v1);
            end
            polyList = [polyList {poly}]; %#ok<AGROW>
        end
    end
end


%% topology management
methods
    function [b1, b2] = isManifold(obj)
        % Check whether the mesh may be considered as manifold.
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
        %     mesh = TopologicalTriMesh3D.createOctahedron;
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
    
    function res = boundary(obj)
        % Compute the boundary of this mesh as a new mesh (can be empty).
        
        edgeInds = boundaryEdgeIndices(obj);
        
        vertexInds = unique(obj.Edges(edgeInds, :));
        
        newVertices = obj.Vertices(vertexInds, :);
        
        newEdges = obj.Edges(edgeInds, :);

        % recompute vertex labels
        for i = 1:numel(newEdges)
            newEdges(i) = find(vertexInds == newEdges(i));
        end
        
        % create the resulting mesh
        res = TopologicalTriMesh3D();
        setVertices(res, newVertices);
        res.Edges = newEdges;
    end
    
    function inds = boundaryEdgeIndices(obj)
        % Find boundary edges and returns their indices.
        
        ensureValidEdges(obj);
        ensureValidEdgeFaces(obj);
   
        % identifies edges adjacent to exactly 1 face.
        inds = find(cellfun(@(x) length(x) == 1, obj.EdgeFaces));
    end
    
    function b = isBoundaryEdge(obj, edgeInd)
        % Check if the  specified edge is boundary.
        
        ensureValidEdgeFaces(obj);
        
        % returns true if edge is adjacent to exactly one face
        b = length(obj.EdgeFaces{edgeInd}) == 1;
    end
    
    function b = isBoundaryVertex(obj, vertexInd)
        % Check if the  specified edge is boundary.
        
        ensureValidEdgeFaces(obj);
        ensureValidVertexEdges(obj);
        
        % returns true if vertex is adjacent to at least one boundary edge
        edgeInds = obj.VertexEdges{vertexInd};
        b = cellfun(@(x) length(x) == 1, obj.EdgeFaces(edgeInds));
        b = any(b);
    end
    
    function res = trimmedMesh(obj)
        % Create new mesh without empty vertices.
        
        % identify vertices referenced by a face
        vertexUsed = false(size(obj.Vertices, 1), 1);
        vertexUsed(unique(obj.Faces(:))) = true;
        inds = find(vertexUsed);
        
        % compute map from old index to new index
        newInds = zeros(size(obj.Vertices, 1), 1);
        for iIndex = 1:length(inds)
            newInds(inds(iIndex)) = iIndex;
        end
        
        % change labels of vertices referenced by faces
        faces2 = zeros(size(obj.Faces));
        for i = 1:numel(faces2)
            faces2(i) = find(inds == obj.Faces(i), 1);
        end
        
        % create new mesh
        res = TopologicalTriMesh3D(obj.Vertices(inds, :), faces2);
    end
end

%% Vertex management methods
methods
    function nv = vertexCount(obj)
        % Return the number of valid vertices.
        % (the result may be different from the size of the Vertices array)
        %
        % NV = vertexCount(MESH);
        %
        % See Also
        %   size, faceCount, edgeCount
        
        nv = sum(obj.VertexValidities);
    end
    
    function ind = addVertex(obj, position)
        if any(size(position) ~= [1 3])
            error('Require a 1-by-3 array of coordinates as input argument');
        end
        obj.Vertices = [obj.Vertices ; position];
        obj.VertexValidities = [obj.VertexValidities ; true];
        ind = size(obj.Vertices, 1);

        % optionnally updates topological data structures
        if ~isempty(obj.VertexEdges)
            obj.VertexEdges{ind} = [];
        end
        if ~isempty(obj.VertexFaces)
            obj.VertexFaces{ind} = [];
        end
    end
    
    function removeVertex(obj, indVertex)
        
        if ~isscalar(indVertex)
            error('vertex index must be a scalar');
        end
        if indVertex > size(obj.Vertices, 1)
            error('vertex index %d is greater than vertex number', indVertex);
        end
        if ~obj.VertexValidities(indVertex)
            error('vertex index %d is already removed', indVertex);
        end
        
        if ~isempty(obj.VertexEdges)
            if ~isempty(obj.VertexEdges{indVertex})
                error('Can not remove vertex %d as it is adjacent to %d edges', ...
                    indVertex, length(obj.VertexEdges{indVertex}));
            end
        end
        if ~isempty(obj.VertexFaces)
            if ~isempty(obj.VertexFaces{indVertex})
                error('Can not remove vertex %d as it is adjacent to %d faces', ...
                    indVertex, length(obj.VertexFaces{indVertex}));
            end
        end
        if any(any(obj.Faces == indVertex))
             error('Can not remove vertex %d as it is adjacent to at least one face', ...
                 indVertex);
        end
        if ~isempty(obj.Edges)
            if any(any(obj.Edges == indVertex))
                error('Can not remove vertex %d as it is adjacent to at least one edge', ...
                    indVertex);
            end
        end
        obj.VertexValidities(indVertex) = false;
    end
    
    function setVertices(obj, coords)
        % Replaces the array of vertices with the new position array.
        obj.Vertices = coords;
        obj.VertexValidities = true(size(obj.Vertices, 1), 1);
    end
end

%% Edge management methods
methods
    function ne = edgeCount(obj)
        % Count the number of edges within the mesh.
        %
        % NE = edgeCount(MESH);
        %
        % See Also
        %   size, faceCount, vertexCount
        
        if isempty(obj.Edges)
            computeEdges(obj);
        end
        ne = sum(sum(obj.Edges == -1, 2) == 0);
    end
    
    function ind = addEdge(obj, edgeVertices)
        % Add an edge given vertex indices and return new edge index.
        
        % identify extremity vertex indices
        edge = sort(edgeVertices, 2); 
        v1 = edge(1); 
        v2 = edge(2);
        
        % check if edge exists
        if ~isempty(obj.Edges)
            if ~isempty(find(obj.Edges(:,1) == v1 & obj.Edges(:,2) == v2, 1))
                error('Edge is already into mesh');
            end
        end
        
        % add the edge
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
    function nf = faceCount(obj)
        % Number of faces within the mesh.
        % Can be different from the size of the Faces array, as some faces
        % may be marked as invalid.
        %
        % Usage
        %   NF = faceCount(MESH);
        %
        % See Also
        %   size, vertexCount, edgeCount

        nf = sum(sum(obj.Faces < 1, 2) == 0);
    end
    
    function indFace = addFace(obj, vertexInds)
        % Add a face given 3 vertex indices and return new face index.
        
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
                obj.EdgeFaces{edgeInds(1)} = [obj.EdgeFaces{edgeInds(1)} indFace];
                obj.EdgeFaces{edgeInds(2)} = [obj.EdgeFaces{edgeInds(2)} indFace];
                obj.EdgeFaces{edgeInds(3)} = [obj.EdgeFaces{edgeInds(3)} indFace];
            end
        end
    end
    
    function removeFace(obj, faceIndex)
        % Remove a face. This resets topological properties.
        obj.Faces(faceIndex, :) = [];
        
        % TODO: update edges instead of clearing
        clearEdges(obj);
    end
    
    function adj = adjacencyMatrix(obj)
        % Compute vertex-to-vertex adjacency matrix.

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
        % Reset all data related to edges as well as VertexFaces.
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
        % Update the property "Edges" from the Faces array.
        
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
        % Compute Face-to-Edge index array, as a Nf-by-3 array.
        % (each face is connected to exactly three edges)
        
        ensureValidEdges(obj);
        
        nFaces = faceCount(obj);
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
        % Compute index of faces adjacent to each edge of a mesh.
        %
        %   Compute index array of faces adjacent to each edge of a mesh.
        %   The result EF has as many rows as the number of edges, and two column.
        %   The first column contains index of faces located on the left of the
        %   corresponding edge, whereas the second column contains index of the
        %   face located on the right. Some indices may be 0 if the mesh is not
        %   'closed'.
        
        ensureValidEdges(obj);
        % indices of faces adjacent to each edge
        obj.EdgeFaces = cell(1, edgeCount(obj));
        
        for iFace = 1:faceCount(obj)
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
                    warning('TopologicalTriMesh3D:IllegalTopology', ...
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
    
    function ensureValidVertexEdges(obj)
        if isempty(obj.VertexEdges)
            computeVertexEdges(obj);
        end
    end
    
    function computeVertexEdges(obj)
        % updates the property "VertexEdges".
        ensureValidEdges(obj);
        
        nVertices = size(obj.Vertices, 1);
        obj.VertexEdges = cell(1, nVertices);
        for iEdge = 1:size(obj.Edges, 1)
            v1 = obj.Edges(iEdge, 1);
            inds = obj.VertexEdges{v1};
            inds = [inds iEdge]; %#ok<AGROW>
            obj.VertexEdges{v1} = inds;
            v2 = obj.Edges(iEdge, 2);
            inds = obj.VertexEdges{v2};
            inds = [inds iEdge]; %#ok<AGROW>
            obj.VertexEdges{v2} = inds;
        end
    end
    
end % end methods


%% Geometric transforms
methods
    function res = scale(obj, varargin)
        % Returns a scaled version of this geometry.
        factor = varargin{1};
        res = TopologicalTriMesh3D(obj.Vertices * factor, obj.Faces);
    end
    
    function res = translate(obj, varargin)
        % Returns a translated version of this geometry.
        shift = varargin{1};
        res = TopologicalTriMesh3D(bsxfun(@plus, obj.Vertices, shift), obj.Faces);
    end
    
end % end methods

%% Overload Matlab functions 
methods
    function dims = size(obj, varargin)
        % Override the size function to return number of mesh elements.
        %
        % S = size(MESH)
        % Returns a vector row with the number of vertices, of edges, and
        % of faces.
        %
        % S = size(MESH, I)
        % Returns the number of elements according to I, that can be 1, 2
        % or 3.
        %
        if isempty(varargin)
            dims = [size(obj.Vertices, 1) size(obj.Edges, 1) size(obj.Faces, 1)];
        else
            nd = varargin{1};
            if nd == 1
                dims = size(obj.Vertices, 1);
            elseif nd == 2
                dims = size(obj.Edges, 1);
            elseif nd == 3
                dims = size(obj.Faces, 1);
            else
                error('Dimension argument must be between 1 and 3');
            end
        end
    end
end

end % end classdef
