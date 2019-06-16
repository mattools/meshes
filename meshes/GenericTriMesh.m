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
end

%% Geometry methods
methods
    function box = boundingBox(obj)
        % Returns the bounding box of this shape.
        mini = min(obj.Vertices);
        maxi = max(obj.Vertices);
        box = Box3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function res = clipVertices(obj, box)
        
        % clip the vertices
        % get bounding box limits
        xmin = box(1); xmax = box(2);
        ymin = box(3); ymax = box(4);
        zmin = box(5); zmax = box(6);
        
        % compute indices of points inside visible area
        xOk = obj.Vertices(:,1) >= xmin & obj.Vertices(:,1) <= xmax;
        yOk = obj.Vertices(:,2) >= ymin & obj.Vertices(:,2) <= ymax;
        zOk = obj.Vertices(:,3) >= zmin & obj.Vertices(:,3) <= zmax;
        
        % select vertices
        inds = find(xOk & yOk & zOk);
        newVertices = obj.Vertices(inds, :);
        
        % create index array for face indices relabeling
        refInds = zeros(1, length(xOk));
        for i = 1:length(inds)
            refInds(inds(i)) = i;
        end
        
        % select the faces with all vertices within the box
        indFaces = sum(~ismember(obj.Faces, inds), 2) == 0;
        newFaces = refInds(obj.Faces(indFaces, :));
        
        res = GenericTriMesh(newVertices, newFaces);
    end
    
    function centroids = faceCentroids(obj)
        % computes the normals of all faces in the mesh
        
        nf = size(obj.Faces, 1);
        centroids = zeros(nf, 3);
        % For triangular meshes, uses accelerated method
        % (taken from https://github.com/alecjacobson/gptoolbox)
        for ff = 1:3
            centroids = centroids + obj.Vertices(obj.Faces(:,ff),:) / 3.0;
        end
    end
    
    function normals = faceNormals(obj)

        % compute vector of first edges
        v1 = obj.Vertices(obj.Faces(:,1),:);
        v12 = obj.Vertices(obj.Faces(:,2),:) - v1;
        v13 = obj.Vertices(obj.Faces(:,3),:) - v1;
        
        % compute normals using cross product (vectors have same size)
        normals = cross(v12, v13, 2);
    end
    
    function normals = vertexNormals(obj)
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
        v1 = obj.Vertices(obj.Edges(edgeIndex,1), :);
        v2 = obj.Vertices(obj.Edges(edgeIndex,2), :);
        centro = (v1 + v2) / 2;
    end
    
    function [dist, proj] = distanceToPoint(obj, point, varargin)
        % Shortest distance between a (3D) point and the mesh.
        %
        %   DIST = distanceToPoint(OBJ, POINT)
        %   Returns the shortest distance between the query point POINT and the
        %   triangular mesh.
        %
        %   [DIST, PROJ] = distanceToPoint(...)
        %   Also returns the projection of the query point on the triangular mesh.
        %
        %
        %   Example
        %     lx = linspace(-1, 1, 100);
        %     [x, y, z] = meshgrid(lx, lx, .5);
        %     pts = [x(:) y(:) z(:)];
        %     dists = ico3s.distanceToPoint(pts);
        %     distMap = reshape(dists, size(x));
        %     figure; imshow(distMap);colormap jet; colorbar
        %
        %   References
        %   * "Distance Between Point and Triangle in 3D", David Eberly (1999)
        %   https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
        %   * <a href="matlab:
        %     web('https://fr.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d')
        %   ">Distance between a point and a triangle in 3d</a>, by Gwendolyn Fischer.
        %   * <a href="matlab:
        %     web('https://fr.mathworks.com/matlabcentral/fileexchange/52882-point2trimesh------distance%C2%A0between-point-and-triangulated-surface')
        %   ">Distance Between Point and Triangulated Surface</a>, by Daniel Frisch.
        
        % Vectorized version of the distanceToPoint function
        %
        %   output = distancePointTrimesh_vectorized(input)
        %
        %   This version is  vectorized over faces: for each query point, the
        %   minimum distance to each triangular face is computed in parallel.
        %   Then the minimum distance over faces is kept.
        %
        
        % Regions are not numbered as in the original paper of D. Eberly to allow
        % automated computation of regions from the 3 conditions on lines.
        % Region indices are computed as follow:
        %   IND = b2 * 2^2 + b1 * 2 + b0
        % with:
        %   b0 = 1 if s < 0, 0 otherwise
        %   b1 = 1 if t < 0, 0 otherwise
        %   b2 = 1 if s+t > 1, 0 otherwise
        % resulting ion the following region indices:
        %        /\ t
        %        |
        %   \ R5 |
        %    \   |
        %     \  |
        %      \ |
        %       \| P3
        %        *
        %        |\
        %        | \
        %   R1   |  \   R4
        %        |   \
        %        | R0 \
        %        |     \
        %        | P1   \ P2
        %  ------*-------*------> s
        %        |        \
        %   R3   |   R2    \   R6
        
        % allocate memory for result
        nPoints = size(point, 1);
        dist = zeros(nPoints, 1);
        proj = zeros(nPoints, 3);
        
        % triangle origins and direction vectors
        p1  = obj.Vertices(obj.Faces(:,1),:);
        v12 = obj.Vertices(obj.Faces(:,2),:) - p1;
        v13 = obj.Vertices(obj.Faces(:,3),:) - p1;
        
        % identify coefficients of second order equation that do not depend on
        % query point
        a = dot(v12, v12, 2);
        b = dot(v12, v13, 2);
        c = dot(v13, v13, 2);
        
        % iterate on query points
        for i = 1:nPoints
            % coefficients of second order equation that depend on query point
            diffP = bsxfun(@minus, p1, point(i, :));
            d = dot(v12, diffP, 2);
            e = dot(v13, diffP, 2);
            
            % compute position of projected point in the plane of the triangle
            det = a .* c - b .* b;
            s   = b .* e - c .* d;
            t   = b .* d - a .* e;
            
            % compute region index (one for each face)
            regIndex = (s < 0) + 2 * (t < 0) + 4 * (s + t > det);
            
            % for each region, process all faces whose projection fall within it
            
            % region 0
            % the minimum distance occurs inside the triangle
            inds = regIndex == 0;
            s(inds) = s(inds) ./ det(inds);
            t(inds) = t(inds) ./ det(inds);
            
            % region 1 (formerly region 3)
            % The minimum distance must occur on the line s = 0
            inds = find(regIndex == 1);
            s(inds) = 0;
            t(inds(e(inds) >= 0)) = 0;
            inds2 = inds(e(inds) < 0);
            bool3 = c(inds2) <= -e(inds2);
            t(inds2(bool3)) = 1;
            inds3 = inds2(~bool3);
            t(inds3) = -e(inds3) ./ c(inds3);
            
            % region 2 (formerly region 5)
            % The minimum distance must occur on the line t = 0
            inds = find(regIndex == 2);
            t(inds) = 0;
            s(inds(d(inds) >= 0)) = 0;
            inds2 = inds(d(inds) < 0);
            bool3 = a(inds2) <= -d(inds2);
            s(inds2(bool3)) = 1;
            inds3 = inds2(~bool3);
            s(inds3) = -d(inds3) ./ a(inds3);
            
            % region 3 (formerly region 4)
            % The minimum distance must occur
            % * on the line t = 0
            % * on the line s = 0 with t >= 0
            % * at the intersection of the two lines
            inds = find(regIndex == 3);
            inds2 = inds(d(inds) < 0);
            % minimum on edge t = 0 with s > 0.
            t(inds2) = 0;
            bool3 = a(inds2) <= -d(inds2);
            s(inds2(bool3)) = 1;
            inds3 = inds2(~bool3);
            s(inds3) = -d(inds3) ./ a(inds3);
            inds2 = inds(d(inds) >= 0);
            % minimum on edge s = 0
            s(inds2) = 0;
            bool3 = e(inds2) >= 0;
            t(inds2(bool3)) = 0;
            bool3 = e(inds2) < 0 & c(inds2) <= e(inds2);
            t(inds2(bool3)) = 1;
            bool3 = e(inds2) < 0 & c(inds2) > e(inds2);
            inds3 = inds2(bool3);
            t(inds3) = -e(inds3) ./ c(inds3);
            
            % region 4 (formerly region 1)
            % The minimum distance must occur on the line s + t = 1
            inds = find(regIndex == 4);
            numer = (c(inds) + e(inds)) - (b(inds) + d(inds));
            s(inds(numer <= 0)) = 0;
            inds2 = inds(numer > 0);
            numer = numer(numer > 0);
            denom = a(inds2) - 2 * b(inds2) + c(inds2);
            s(inds2(numer > denom)) = 1;
            bool3 = numer <= denom;
            s(inds2(bool3)) = numer(bool3) ./ denom(bool3);
            t(inds) = 1 - s(inds);
            
            % Region 5 (formerly region 2)
            % The minimum distance must occur:
            % * on the line s + t = 1
            % * on the line s = 0 with t <= 1
            % * or at the intersection of the two (s=0; t=1)
            inds = find(regIndex == 5);
            tmp0 = b(inds) + d(inds);
            tmp1 = c(inds) + e(inds);
            % minimum on edge s+t = 1, with s > 0
            bool2 = tmp1 > tmp0;
            inds2 = inds(bool2);
            numer = tmp1(bool2) - tmp0(bool2);
            denom = a(inds2) - 2 * b(inds2) + c(inds2);
            bool3 = numer < denom;
            s(inds2(~bool3)) = 1;
            inds3 = inds2(bool3);
            s(inds3) = numer(bool3) ./ denom(bool3);
            t(inds2) = 1 - s(inds2);
            % minimum on edge s = 0, with t <= 1
            inds2 = inds(~bool2);
            s(inds2) = 0;
            t(inds2(tmp1(~bool2) <= 0)) = 1;
            t(inds2(tmp1(~bool2) > 0 & e(inds2) >= 0)) = 0;
            inds3 = inds2(tmp1(~bool2) > 0 & e(inds2) < 0);
            t(inds3) = -e(inds3) ./ c(inds3);
            
            % region 6 (formerly region 6)
            % The minimum distance must occur
            % * on the line s + t = 1
            % * on the line t = 0, with s <= 1
            % * at the intersection of the two lines (s=1; t=0)
            inds = find(regIndex == 6);
            tmp0 = b(inds) + e(inds);
            tmp1 = a(inds) + d(inds);
            % minimum on edge s+t=1, with t > 0
            bool2 = tmp1 > tmp0;
            inds2 = inds(bool2);
            numer = tmp1(bool2) - tmp0(bool2);
            denom = a(inds2) - 2 * b(inds2) + c(inds2);
            bool3 = numer <= denom;
            t(inds2(~bool3)) = 1;
            inds3 = inds2(bool3);
            t(inds3) = numer(bool3) ./ denom(bool3);
            s(inds2) = 1 - t(inds2);
            % minimum on edge t = 0 with s <= 1
            inds2 = inds(~bool2);
            t(inds2) = 0;
            s(inds2(tmp1(~bool2) <= 0)) = 1;
            s(inds2(tmp1(~bool2) > 0 & d(inds2) >= 0)) = 0;
            inds3 = inds2(tmp1(~bool2) > 0 & d(inds2) < 0);
            s(inds3) = -d(inds3) ./ a(inds3);
            
            % compute coordinates of closest point on plane
            projList = p1 + bsxfun(@times, s, v12) + bsxfun(@times, t, v13);
            
            % squared distance between point and closest point on plane
            [dist(i), ind] = min(sum((bsxfun(@minus, point(i,:), projList)).^2, 2));
            
            % keep the valid projection
            proj(i, :) = projList(ind,:);
        end
        
        % convert squared distance to distance
        dist = sqrt(dist);
    end
end


%% Topological queries
methods
    function res = vertexLink(obj, vertexIndex)
        % returns the link around the specifid vertex as a new mesh
        
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
        res = GenericTriMesh();
        res.Vertices = newVertices;
        res.Edges = newEdges;
    end
    
    function polyList = vertexLinkPolygons(obj, vertexIndex)
        % returns the link around a vertex as a list of 3D polylines
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
    
    function res = boundary(obj)
        % boundary of this mesh as a new mesh (can be empty)
        
        edgeInds = boundaryEdgeIndices(obj);
        
        vertexInds = unique(obj.Edges(edgeInds, :));
        
        newVertices = obj.Vertices(vertexInds, :);
        
        newEdges = obj.Edges(edgeInds, :);

        % recompute vertex labels
        for i = 1:numel(newEdges)
            newEdges(i) = find(vertexInds == newEdges(i));
        end
        
        % create the resulting mesh
        res = GenericTriMesh();
        res.Vertices = newVertices;
        res.Edges = newEdges;
    end
    
    function inds = boundaryEdgeIndices(obj)
        % finds boundary edges and returns their indices
        
        ensureValidEdges(obj);
        ensureValidEdgeFaces(obj);
   
        % identifies edges adjacent to exactly 1 face.
        inds = find(cellfun(@(x) length(x) == 1, obj.EdgeFaces));
    end
    
    function b = isBoundaryEdge(obj, edgeInd)
        % checks if the  specified edge is boundary
        
        ensureValidEdgeFaces(obj);
        
        % returns true if edge is adjacent to exactly one face
        b = length(obj.EdgeFaces{edgeInd}) == 1;
    end
    
    function b = isBoundaryVertex(obj, vertexInd)
        % checks if the  specified edge is boundary
        
        ensureValidEdgeFaces(obj);
        ensureValidVertexEdges(obj);
        
        % returns true if vertex is adjacent to at least one boundary edge
        edgeInds = obj.VertexEdges{vertexInd};
        b = cellfun(@(x) length(x) == 1, obj.EdgeFaces(edgeInds));
        b = any(b);
    end
    
    function res = trimmedMesh(obj)
        % new mesh without empty vertices
        
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
        res = GenericTriMesh(obj.Vertices(inds, :), faces2);
    end
    
    function ccList = connectedComponents(obj)
        % set of edge-connected components
        
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
            cc = trimmedMesh(GenericTriMesh(obj.Vertices, newFaces));
            ccList = [ccList {cc}]; %#ok<AGROW>
        end
        
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
                obj.EdgeFaces{edgeInds(1)} = [obj.EdgeFaces{edgeInds(1)} indFace];
                obj.EdgeFaces{edgeInds(2)} = [obj.EdgeFaces{edgeInds(2)} indFace];
                obj.EdgeFaces{edgeInds(3)} = [obj.EdgeFaces{edgeInds(3)} indFace];
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
        
        ensureValidEdges(obj);
        
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
        
        ensureValidEdges(obj);
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
    
    function ensureValidVertexEdges(obj)
        if isempty(obj.VertexEdges)
            computeVertexEdges(obj);
        end
    end
    
    function computeVertexEdges(obj)
        % updates the property "VertexEdges"
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

%% Drawing functions
methods
    function h = draw(varargin)
        % Draw the faces of this mesh, using the patch function.
        % see also
        %   drawEdges
        
        hh = drawFaces(varargin{:});
                
        if nargout > 0
            h = hh;
        end
    end
    
    function h = drawVertices(varargin)
        % Draw the vertices of this complex
        
        % extract handle of axis to draw in
        if numel(varargin{1}) == 1 && ishghandle(varargin{1}, 'axes')
            hAx = varargin{1};
            varargin(1) = [];
        else
            hAx = gca;
        end

        % extract the mesh instance from the list of input arguments
        obj = varargin{1};
        varargin(1) = [];
        
        % add default drawing options
        inds = [];
        options = {'linestyle', 'none', 'marker', 'o'};

        % extract optional drawing options
        if nargin > 1 && isnumeric(varargin{1})
            % get index of edges to draw
            inds = varargin{1};
            varargin(1) = [];
        end
        if ~isempty(varargin) && ischar(varargin{1})
            if length(varargin) > 1
                options = [options varargin];
            else
                options = varargin;
            end
        end
        
        % Draw 3D points
        if isempty(inds)
            x = obj.Vertices(:, 1);
            y = obj.Vertices(:, 2);
            z = obj.Vertices(:, 3);
        else
            x = obj.Vertices(inds, 1);
            y = obj.Vertices(inds, 2);
            z = obj.Vertices(inds, 3);
        end
        hh = plot3(hAx, x, y, z, options{:});

        % optionnally add style processing
        if ~isempty(varargin) && isa(varargin{1}, 'Style')
            apply(varargin{1}, hh);
        end
                
        if nargout > 0
            h = hh;
        end
    end

    function h = drawEdges(varargin)
        % Draw the edges of this complex
        
        % extract handle of axis to draw in
        if numel(varargin{1}) == 1 && ishghandle(varargin{1}, 'axes')
            hAx = varargin{1};
            varargin(1) = [];
        else
            hAx = gca;
        end

        % extract the mesh instance from the list of input arguments
        obj = varargin{1};
        varargin(1) = [];
        
        % add default drawing options
        inds = [];
        options = {'Color', [0 0 0]};

        % extract optional drawing options
        if nargin > 1 && isnumeric(varargin{1})
            % get index of edges to draw
            inds = varargin{1};
            varargin(1) = [];
        end
        if nargin > 1 && ischar(varargin{1})
            if nargin > 2
                options = [options varargin];
            else
                options = varargin;
            end
        end
        
        % Draw 3D edges
        if isempty(inds)
            x = [obj.Vertices(obj.Edges(:,1), 1) obj.Vertices(obj.Edges(:,2), 1)]';
            y = [obj.Vertices(obj.Edges(:,1), 2) obj.Vertices(obj.Edges(:,2), 2)]';
            z = [obj.Vertices(obj.Edges(:,1), 3) obj.Vertices(obj.Edges(:,2), 3)]';
        else
            x = [obj.Vertices(obj.Edges(inds,1), 1) obj.Vertices(obj.Edges(inds,2), 1)]';
            y = [obj.Vertices(obj.Edges(inds,1), 2) obj.Vertices(obj.Edges(inds,2), 2)]';
            z = [obj.Vertices(obj.Edges(inds,1), 3) obj.Vertices(obj.Edges(inds,2), 3)]';
        end
        hh = plot3(hAx, x, y, z, options{:});

        % optionnally add style processing
        if ~isempty(varargin) && isa(varargin{1}, 'Style')
            apply(varargin{1}, hh);
        end
                
        if nargout > 0
            h = hh;
        end
    end

    function h = drawFaces(varargin)
        % Draw the faces of this mesh, using the patch function.
        %
        % see also
        %   draw, drawVertices, drawEdges
        
        % extract handle of axis to draw in
        if numel(varargin{1}) == 1 && ishghandle(varargin{1}, 'axes')
            hAx = varargin{1};
            varargin(1) = [];
        else
            hAx = gca;
        end

        % extract the mesh instance from the list of input arguments
        obj = varargin{1};
        varargin(1) = [];
        
        % add default drawing options
        inds = [];
        options = {'FaceColor', [.75 .75 .75]};

        % extract optional drawing options
        if nargin > 1 && isnumeric(varargin{1})
            % get index of faces to draw
            inds = varargin{1};
            varargin(1) = [];
        end
        if nargin > 1 && ischar(varargin{1})
            options = [options varargin];
        end
        
        if length(options) == 1
            options = [{'facecolor', [.75 .75 .75]} options];
        end

        if isempty(inds)
            hh = patch('Parent', hAx, ...
                'vertices', obj.Vertices, 'faces', obj.Faces, ...
                options{:} );
        else
            hh = patch('Parent', hAx, ...
                'vertices', obj.Vertices, 'faces', obj.Faces(inds, :), ...
                options{:} );
        end
        
        % optionnally add style processing
        if ~isempty(varargin) && isa(varargin{1}, 'Style')
            apply(varargin{1}, hh);
        end
                
        if nargout > 0
            h = hh;
        end
    end
    
    function h = drawFaceNormals(obj, varargin)
        
        % compute vector data
        c = faceCentroids(obj);
        n = faceNormals(obj);
        
        % display an arrow for each normal
        hq = quiver3(c(:,1), c(:,2), c(:,3), n(:,1), n(:,2), n(:,3));
        
        % format output
        if nargout > 0
            h = hq;
        end
    end
    
end

%% Geometric transforms
methods
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

%% Overload Matlab functions 
methods
    function dims = size(obj, varargin)
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

