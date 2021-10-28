function h = drawFaces(varargin)
% Draw the faces of this mesh, using the patch function.
%
% see also
%   draw, drawVertices, drawEdges

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

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