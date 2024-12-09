% compute edges for integral torque, splits edges into different
% directions

function [edges,thetaI] = computeEdge3( img, edgeS, num_orient )
% computeEdge3 - Compute and splits edges into 8 different directions for torque, 
% v3 masks the final edges using edgeS and allows one to set the number of orientation
%
% USAGE        
%   edges = computeEdge3( img, edgeS, num_orient )
%
% INPUTS
%   img          - [HxWxD] input RGB or grayscale image
%   edgeS        - [HxW] input edge map
%   num_orient   - (8) number of orientations desired
%
% OUTPUTS
%   edges  	- [HxWxnum_orient] edges detected at 8 orientations
%   thetaI  - Orientations of edges
%
% Copyright (C) 2015 Ching L Teo, University of Maryland College Park, [cteo-at-cs.umd.edu]
%
%    You can redistribute and/or modify this software for non-commercial use
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    For commercial use, contact the author for licensing options. 
%
% Please email me if you find bugs, or have suggestions or questions!

if nargin < 3
    num_orient=8;
end
   
[~, ~, nDim] = size(img);

if nDim == 3
    imGray = rgb2gray( img );
else
    imGray = img;
end


[gx,gy] = gradient( imGray );
    
theta = atan2( gx, gy );
theta( gx==0 & gy==0 ) = NaN;
thetaI=theta;
thetaI(~edgeS) = NaN;

% split edges into num_orient directions
theta = theta + pi/num_orient;
theta(theta<0) = theta(theta<0) + 2.0*pi;
edges = zeros([size(edgeS),num_orient]);
for i=1:num_orient
    mask = ( theta>=(i-1)/num_orient*2.0*pi & theta<i/num_orient*2.0*pi );
    edges(:,:,i) = edgeS .* mask;
end

