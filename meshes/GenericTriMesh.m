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


%% Local methods for updating local Properties 
methods
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
        res = TriMesh3D(obj.VertexCoords * factor, obj.FaceVertexInds);
    end
    
    function res = translate(obj, varargin)
        % Returns a translated version of this geometry
        shift = varargin{1};
        res = TriMesh3D(bsxfun(@plus, obj.vertexCords, shift), obj.FaceVertexInds);
    end
    
end % end methods


%% Global use methods
methods
    function disp(obj)
        nv = size(obj.Vertices, 1);
        ne = size(obj.Edges, 1);
        nf = size(obj.Faces, 1);
        fprintf('Generic Mesh with %d vertices, %d edges, %d faces\n', ...
            nv, ne, nf);
    end
    
end % end methods



end % end classdef

