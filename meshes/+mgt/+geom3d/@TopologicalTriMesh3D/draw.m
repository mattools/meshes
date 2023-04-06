function h = draw(varargin)
% Draw the faces of this mesh, using the patch function.
%
%   output = truc(input)
%
%   Example
%   truc
%
%   See also
%     drawFaces, drawEdges
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-10-28,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE.

hh = drawFaces(varargin{:});

if nargout > 0
    h = hh;
end
