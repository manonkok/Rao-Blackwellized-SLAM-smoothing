function h = imagescalpha(x,y,C,A,alim)
% IMAGESCALPHA - Displays an image with scaled alpha values
%
% Syntax:
%   h = imagescalpha(x, y, C, A, alim)
%
% In:
%   x    - X-coordinates for the image
%   y    - Y-coordinates for the image
%   C    - Color data for the image
%   A    - Alpha (transparency) data for the image
%   alim - Limits for scaling the alpha data (optional)
%
% Out:
%   h    - Handle to the image object
%
% Description:
%   This function displays an image with its alpha (transparency) values 
%   scaled according to the provided limits. The function scales the alpha 
%   values within the specified range, allowing for adjustable transparency 
%   effects. If 'alim' is not provided, the function automatically 
%   determines the limits based on the alpha data.
%
% See also:
%   imagesc
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

% Set automatic limits
if nargin < 5 || isempty(alim)
  alim = [min(A(:)) max(A(:))];
end

%% Scale alpha

A = A-alim(1);
A = 1 - A/(alim(2)-alim(1));
A = uint8(255*A);

%% Visualize

h=imagesc(x,y,C);
set(h,'AlphaData',A)
set(h,'AlphaDataMapping','direct')











