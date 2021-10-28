function h = drawEdges(varargin)
% Draw the edges of this mesh.
%
%   output = truc(input)
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
