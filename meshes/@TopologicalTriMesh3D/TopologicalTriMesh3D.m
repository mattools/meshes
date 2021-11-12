classdef TopologicalTriMesh3D < handle
% A 3D triangular mesh that stores the full topology information.
%
%   The TopologicalTriMesh3D class is an implementation of 3D triangular
%   mesh that stores the full topology information between vertices, edges
%   and faces. It should be possible to retrieve with edges are incident to
%   a vertex, which faces are incident to an edge, and so on...
%   It will be possible to remove elements from the mesh.
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
    % The array of vertex coordinates, as a Nv-by-3 numeric array.
    Vertices;
    
    % The list of edges adjacent to each vertex, as a 1-by-NV cell array.
    VertexEdges = cell(1, 0);
    
    % The list of faces adjacent to each vertex, as a 1-by-NV cell array.
    VertexFaces = cell(1, 0);
    
    % An array of boolean flags indicating for which vertex if it is valid.
    % (default is true for each vertex)
    ValidVertices = true([0 1]);
    
    
    % The topology on edges, as a NE-by-2 array containing vertex indices.
    Edges = zeros([0 2]);
    
    % Cell array containing indices of faces associated to each edge.
    % should be 2 for regular edges, or 1 for boundary edges, but can be
    % more for non manifold edges.
    EdgeFaces = cell(1, 0);

    % An array of boolean flags indicating for which edge if it is valid.
    % (default is true for each edge)
    ValidEdges = true([0 1]);
    
    
    % The NF-by-3 array containing vertex indices of each face.
    % Some rows may contains [-1 -1 -1] if the face was removed.
    Faces = zeros([0 3]);
    
    % The array of edges adjacent to each face, as a Nf-by-3 array of edge
    % indices.
    FaceEdges = zeros([0 3]);
    
    % An array of boolean flags indicating for which face if it is valid.
    % (default is true for each face)
    ValidFaces = true([0 1]);
    
end % end properties

%% Static factories
methods (Static)
    function obj = createOctahedron()
        vertices = [1 0 0;0 1 0;-1 0 0;0 -1 0;0 0 1;0 0 -1];
        faces = [1 2 5;2 3 5;3 4 5;4 1 5;1 6 2;2 6 3;3 6 4;1 4 6];
        obj = TopologicalTriMesh3D(vertices, faces);
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
            initializeVertices(obj, var1.Vertices);
            obj.Faces = var1.Faces;
            
        elseif nargin == 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
            % classical vertex-face constructor
            initializeVertices(obj, varargin{1});
            
            faces = varargin{2};
            for iFace = 1:size(faces, 1)
                addFace(obj, faces(iFace, :));
            end
            
        else
            error('Unable to interpret input arguments');
            
        end
        
        function initializeVertices(obj, vertices)
            % Initialize vertices and vertex-related data.
            nv = size(vertices, 1);
            obj.Vertices = vertices;
            obj.ValidVertices = true(nv, 1);
            obj.VertexEdges = cell(1, nv);
            obj.VertexFaces = cell(1, nv);
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
        
        % TODO: updates faceEdges info
        % TODO: updates vertexEdges info
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
        
        % TODO: updates Edges info
        % TODO: updates faceEdges info
        % TODO: updates vertexEdges info
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
        
        % TODO: updates faceEdges info
        % TODO: updates vertexEdges info
        
    end
end

%% Geometry methods
methods
    function box = bounds(obj)
        % Return the bounds of this mesh.
        vertices = obj.Vertices(obj.ValidVertices, :);
        mini = min(vertices);
        maxi = max(vertices);
        box = Bounds3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function [lengths, inds] = edgeLength(obj, varargin)
        % Compute length of edges in mesh.
        %
        % LEN = edgeLength(MESH)
        % LEN = edgeLength(MESH, EINDS)
        % [LEN, INDS] = edgeLength(...)
        % also returns the indices of edges for which length was computed.
        
        if isempty(varargin)
            inds = find(obj.ValidEdges);
        else
            inds = varargin{1};
        end
        iv1 = obj.Edges(inds,1);
        iv2 = obj.Edges(inds,2);
        lengths = sqrt(sum((obj.Vertices(iv1, :) - obj.Vertices(iv2, :)).^2, 2));
    end
    
    function centroids = faceCentroid(obj, varargin)
        % Compute the normals of all faces in the mesh
        
        if nargin > 1
            % process specified faces
            inds = varargin{1};
        else
            % process only valid faces
            inds = find(obj.ValidFaces);
        end
        
        nf = length(inds);
        centroids = zeros(nf, 3);
        % For triangular meshes, uses accelerated method
        % (taken from https://github.com/alecjacobson/gptoolbox)
        for ff = 1:3
            centroids = centroids + obj.Vertices(obj.Faces(inds,ff),:) / 3.0;
        end
    end
    
    function normals = faceNormal(obj)
        % Compute the normals of faces.

        % process only valid faces
        inds = obj.ValidFaces;
        
        % compute vector of first edges
        v1 = obj.Vertices(obj.Faces(inds,1),:);
        v12 = obj.Vertices(obj.Faces(inds,2),:) - v1;
        v13 = obj.Vertices(obj.Faces(inds,3),:) - v1;
        
        % compute normals using cross product (vectors have same size)
        normals = cross(v12, v13, 2);
        %TODO: return Vector3D?
    end
    
    function normals = vertexNormal(obj)
        % Compute the normals at vertices.
        nv = vertexCount(obj);
        nf = faceCount(obj);
        
        % unit normals to the faces
        fNormals = faceNormals(obj);
        
        % compute normal of each vertex: sum of normals to each face
        normals = zeros(nv, 3);
        inds = find(obj.ValidFaces);
        for i = 1:nf
            face = obj.Faces(inds(i), :);
            for j = 1:length(face)
                v = face(j);
                normals(v, :) = normals(v,:) + fNormals(i,:);
            end
        end
    end

    function centro = edgeCentroid(obj, varargin)
        % Compute centroid of edge.
        
        if isempty(varargin)
            inds = obj.ValidEdges;
        else
            inds = varargin{1};
        end
        v1 = obj.Vertices(obj.Edges(inds,1), :);
        v2 = obj.Vertices(obj.Edges(inds,2), :);
        centro = (v1 + v2) / 2;
    end
end


%% Topological queries
methods
    function res = vertexLink(obj, vertexIndex)
        % Return the 1-link around the specifid vertex as a new mesh.
        
        % get index of edges adjacent to the vertex
        edgeInds = obj.VertexEdges{vertexIndex};
        % index of adjacent vertices
        adjVertexInds = unique(obj.Edges(edgeInds, :));
        adjVertexInds(adjVertexInds == vertexIndex) = [];
        
        % find index of edges belonging to the link.
        % need to iterate over adjacent vertices, and identify edges
        % containing two adjcent vertices
        newEdgeInds = [];
        for iv = 1:length(adjVertexInds)
            vertexIdx = adjVertexInds(iv);
            % index of edges incident to current vertex
            edgeInds = obj.VertexEdges{vertexIdx};
            % vertex edges indices of edges around current vertex
            vertexEdgeArray = obj.Edges(edgeInds, :);
            inds2 = sum(ismember(vertexEdgeArray, adjVertexInds), 2) == 2;
            newEdgeInds = [newEdgeInds ; edgeInds(inds2)']; %#ok<AGROW>
        end
        newEdgeInds = unique(newEdgeInds);
        
        % create the resulting mesh
        res = TopologicalTriMesh3D();
        % add each vertex, keeping the mapping of old vertex indices to new
        % vertex indices 
        vertexIndsMap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
        for iv = 1:length(adjVertexInds)
            oldIndex = adjVertexInds(iv);
            newIndex = addVertex(res, obj.Vertices(oldIndex, :));
            vertexIndsMap(oldIndex) = newIndex;
        end
        % add edges
        for ie = 1:length(newEdgeInds)
            edgeIdx = newEdgeInds(ie);
            iv1 = obj.Edges(edgeIdx, 1);
            iv2 = obj.Edges(edgeIdx, 2);
            iv1n = vertexIndsMap(iv1);
            iv2n = vertexIndsMap(iv2);
            addEdge(res, [iv1n iv2n]);
        end
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

        % compute degree of each edge
        nEdges = edgeCount(obj);
        edgeFaceDegrees = zeros(nEdges, 1);
        for iEdge = 1:nEdges
            edgeFaceDegrees(iEdge) = length(obj.EdgeFaces{iEdge});
        end
        
        % for each face, concatenate the face degree of each edge
        nFaces = faceCount(obj);
        faceEdgeDegrees = zeros(nFaces, 3);
        for iFace = 1:nFaces
            edgeInds = obj.FaceEdges(iFace, :);
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
        
        % identifies edges adjacent to exactly 1 face.
        inds = find(cellfun(@(x) length(x) == 1, obj.EdgeFaces));
    end
    
    function b = isBoundaryEdge(obj, edgeInd)
        % Check if the  specified edge is boundary.
        
        % returns true if edge is adjacent to exactly one face
        b = length(obj.EdgeFaces{edgeInd}) == 1;
    end
    
    function b = isBoundaryVertex(obj, vertexInd)
        % Check if the  specified edge is boundary.
        
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
        
        nv = sum(obj.ValidVertices);
    end
    
    function ind = addVertex(obj, position)
        % Add one or several vertex(ices) from position(s).
        %
        %  addVertex(MESH, [VX VY VZ]);
        %  VIDX = addVertex(MESH, [VX VY VZ]);
        
        % check input validity
        if ~isnumeric(position) || size(position, 2) ~= 3
            error('Position must be a numeric array with three columns');
        end
        
        % previous max vertex index
        ind0 = size(obj.Vertices, 1);
        
        % add the new vertex
        nv = size(position, 1);
        obj.Vertices = [obj.Vertices ; position];
        obj.ValidVertices = [obj.ValidVertices ; true(nv,1)];
        ind = ((ind0+1):size(obj.Vertices, 1))';

        % initialize vertex data
        obj.VertexEdges(ind) = cell(nv, 1);
        obj.VertexFaces(ind) = cell(nv, 1);
    end
    
    function vertexInds = removeInvalidVertices(obj)
        % Remove vertices marked as invalid from the mesh.
        
        % find index of vertices to remove (if any)
        vertexInds = find(~obj.ValidVertices);
        
        % Compute the mapping from old indices to new indices
        vertexIndexMap = zeros(size(obj.Vertices, 1), 1);
        vertexIndexMap(obj.ValidVertices) = 1:sum(obj.ValidVertices);
        
        % re-label edge info
        edges = obj.Edges(obj.ValidEdges, :);
        edges = vertexIndexMap(edges);
        obj.Edges(obj.ValidEdges, :) = edges;
        
        % re-label face info
        faces = obj.Faces(obj.ValidFaces, :);
        faces = vertexIndexMap(faces);
        obj.Faces(obj.ValidFaces, :) = faces;
        
        % remove vertex data
        obj.Vertices(vertexInds, :) = [];
        obj.VertexEdges(vertexInds) = [];
        obj.VertexFaces(vertexInds) = [];
        obj.ValidVertices(vertexInds) = [];
    end
    
    function removeVertex(obj, indVertex)
        
        if ~isscalar(indVertex)
            error('vertex index must be a scalar');
        end
        if indVertex > size(obj.Vertices, 1)
            error('vertex index %d is greater than vertex number', indVertex);
        end
        if ~obj.ValidVertices(indVertex)
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
        obj.ValidVertices(indVertex) = false;
    end
    
    function setVertices(obj, coords)
        % Replaces the array of vertices with the new position array.
        obj.Vertices = coords;
        obj.ValidVertices = true(size(obj.Vertices, 1), 1);
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
        
        ne = sum(obj.ValidEdges);
    end
    
    function ind = addEdge(obj, edgeVertices)
        % Add an edge given vertex indices and return new edge index.
        
        % check input validity
        if any(~obj.ValidVertices(edgeVertices))
            error('Can not create an edge with invalid vertices');
        end
        
        % identify extremity vertex indices
        edge = sort(edgeVertices, 2); 
        v1 = edge(1); 
        v2 = edge(2);

%         % check if edge exists
%         if ~isempty(obj.Edges)
%             if ~isempty(find(obj.Edges(:,1) == v1 & obj.Edges(:,2) == v2, 1))
%                 error('Edge is already into mesh');
%             end
%         end
        
        % add the edge
        obj.Edges = [obj.Edges; edge];
        obj.ValidEdges = [obj.ValidEdges ; true];
        ind = size(obj.Edges, 1);
        
        % add edge index to vertices
        obj.VertexEdges{v1} = [obj.VertexEdges{v1} ind];
        obj.VertexEdges{v2} = [obj.VertexEdges{v2} ind];
        
        % create empty list of faces for this edge
        obj.EdgeFaces{ind} = [];
    end
    
    function ind = findEdgeIndex(obj, edgeVertices)
        % Find the index of an edge from index of its vertices.
        % Return 0 if there is no such edge.
        v1 = edgeVertices(1);
        v2 = edgeVertices(2);
        edgesV1 = obj.VertexEdges{v1};
        for ie = 1:length(edgesV1)
            edgeIdx = edgesV1(ie);
            if ismember(v2, obj.Edges(edgeIdx, :))
                ind = edgeIdx;
                return;
            end
        end
        warning('Could not find edge index between vertices %d and %d', v1, v2);
        ind = 0;
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

        nf = sum(obj.ValidFaces);
    end
    
    function indFace = addFace(obj, vertexInds)
        % Add a face given 3 vertex indices and return new face index.
        %
        % FIDX = addFace(mesh, [V1 V2 V3]);
        %
        
        if any(size(vertexInds) ~= [1 3])
            error('Require an array of three indices as input argument');
        end
        if any(~obj.ValidVertices(vertexInds))
            error('Can not create a face with invalid vertices');
        end
        
        % add the new face, and keep its index
        obj.Faces = [obj.Faces ; vertexInds];
        obj.ValidFaces = [obj.ValidFaces ; true];
        indFace = size(obj.Faces, 1);

        % for each vertex, add a reference to the face
        for iv = 1:3
            indV = vertexInds(iv);
            obj.VertexFaces{indV} = [obj.VertexFaces{indV} indFace];
        end
        
        % retrieve the list of edges adjacent to face vertices
        vertexEdgeInds = unique([obj.VertexEdges{vertexInds}]);
        
        % vertex indices of face edges
        newEdges = sort(vertexInds([1 2; 2 3;1 3]), 2);
        edgeInds = zeros(1,3);
        
        % for each edge, retrieve its index, or create a new edge.
        for iEdge = 1:3
            v1 = newEdges(iEdge, 1);
            v2 = newEdges(iEdge, 2);
            edgeIdx = find(obj.Edges(vertexEdgeInds,1) == v1 & obj.Edges(vertexEdgeInds,2) == v2);
            if isempty(edgeIdx)
                edgeIdx = addEdge(obj, [v1 v2]);
            else
                edgeIdx = vertexEdgeInds(edgeIdx);
            end
            edgeInds(iEdge) = edgeIdx;
        end

        obj.FaceEdges(indFace, :) = edgeInds;
        obj.EdgeFaces{edgeInds(1)} = [obj.EdgeFaces{edgeInds(1)} indFace];
        obj.EdgeFaces{edgeInds(2)} = [obj.EdgeFaces{edgeInds(2)} indFace];
        obj.EdgeFaces{edgeInds(3)} = [obj.EdgeFaces{edgeInds(3)} indFace];
    end
    
    function removeFace(obj, faceIndex)
        % Remove a face. 
        obj.ValidFaces(faceIndex) = false;
        
        % detach face from adjacent edges
        edgeInds = obj.FaceEdges(faceIndex, :);
        for ie = 1:3
            inds = obj.EdgeFaces{edgeInds(ie)};
            inds(inds == faceIndex) = [];
            obj.EdgeFaces{edgeInds(ie)} = inds;
        end
        
        % detach face from adjacent vertices
        vertInds = obj.Faces(faceIndex, :);
        for iv = 1:3
            inds = obj.VertexFaces{vertInds(iv)};
            inds(inds == faceIndex) = [];
            obj.VertexFaces{vertInds(iv)} = inds;
        end
    end
    
    function invalidFaceInds = removeInvalidFaces(obj)
        % Remove faces marked as invalid from the mesh.
        
        % find index of faces to remove (if any)
        invalidFaceInds = find(~obj.ValidFaces);
        validFaceInds = find(obj.ValidFaces);
        
        % Compute the mapping from old indices to new indices
        faceIndexMap = zeros(size(obj.Faces, 1), 1);
        faceIndexMap(validFaceInds) = 1:length(validFaceInds);
        
        % cleanup edges
        validEdgeInds = find(obj.ValidEdges);
        for iEdge = 1:length(validEdgeInds)
            idx = validEdgeInds(iEdge);
            obj.EdgeFaces{idx} = faceIndexMap(obj.EdgeFaces{idx});
        end
        
        % cleanup vertices
        validVertexInds = find(obj.ValidVertices);
        for iVertex = 1:length(validVertexInds)
            idx = validVertexInds(iVertex);
            obj.VertexFaces{idx} = faceIndexMap(obj.VertexFaces{idx});
        end
        
        % keep valid face data
        obj.Faces = obj.Faces(validFaceInds, :);
        obj.FaceEdges = obj.FaceEdges(validFaceInds, :);
        obj.ValidFaces = obj.ValidFaces(validFaceInds);
    end
    
    function adj = adjacencyMatrix(obj)
        % Compute vertex-to-vertex adjacency matrix.

        % forces faces to be floating point array, for sparse function
        faces = obj.Faces(obj.ValidFaces, :);
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
            dims = [vertexCount(obj) edgeCount(obj) faceCount(obj)];
        else
            nd = varargin{1};
            if nd == 1
                dims = vertexCount(obj);
            elseif nd == 2
                dims = edgeCount(obj);
            elseif nd == 3
                dims = faceCount(obj);
            else
                error('Dimension argument must be between 1 and 3');
            end
        end
    end
end

end % end classdef
