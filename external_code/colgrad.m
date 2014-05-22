function Igrad = colgrad(I)
% Gradient image of a multi-spectral image.
%
% It calculates the maximum rate of change in the pixel value (gradient) 
% along the spectral dimensions. It is the implementation of the gradient
% calculation algorithm described by Di Zenzo (1986).
%
% It uses sobel operators to estimate the rate of change along the 
% x and y axis for each channel. eg. RGB image and r,g,b unitary vectors.
% 
% u = r.dR/dx + g.dG/dx + b.dB/dx
% v = r.dR/dy + g.dG/dy + b.dB/dy
% 
% Tensor coefficients:
% gxx = <u,u> = u'.u
% gyy = <v,v> = v'.v
% gxy = <u,v> = u'.v
%
% Input image, I
% --------------
% The input image I should be a 3-D multidimensional array with the third
% dimension corresponding to the spectral components. It works with any
% multispectral image and not just with threebands (eg. RGB).
% The class can be uint8, uint16, int16, single, or double, and it
% must be nonsparse.
% 
% Output image, Igrad
% -------------------
% The output image is the gradient image and is of class double
% 
% Example 1
% ---------
%   I = imread('peppers.png');
%   Igrad = colgrad(I);
%   figure, imshow(Igrad), title('Image gradient: peppers')
%   level = graythresh(Igrad);
%   figure, imshow(Igrad>level), title('Threshold image: peppers')
%
% Example 2
% ---------
%   I = multibandread('paris.lan', [512, 512, 7], 'uint8=>uint8',128,...
%   'bil', 'ieee-le', {'Band','Direct',[1:7]});
%   Igrad = colgrad(I);
%   figure, imshow(Igrad), title('Image gradient: paris')
%   level = graythresh(Igrad);
%   figure, imshow(Igrad>level), title('Threshold image: paris')
%
%   Carlos Gias
%   Date: 21/08/2011

% Reference:
% Di Zenzo, A. Note on the Gradient of a Multi-Image (1986)
% Computer Vision, Graphics, and Image Processing 33, 116-125.

% Normalize the original image between 0-1
max_int = double(max(I(:)));
In = double(I)/max_int;

% Image dimensions
[n_y, n_x, n_ch] = size(In);

% Convolve the image accross channels for horizontal and vertical axis 
% using sobel operators to calculate the partial derivatives
h_hor = 0.25*fspecial('sobel');
h_ver = h_hor';

I_h = zeros(n_y,n_x,n_ch);
I_v = zeros(n_y,n_x,n_ch);
for i=1:n_ch
    I_h(:,:,i) = imfilter(In(:,:,i), h_hor, 'replicate');
    I_v(:,:,i) = imfilter(In(:,:,i), h_ver, 'replicate');
end

% Calculate the tensor coefficients
gxx = sum(I_h.^2, 3);
gyy = sum(I_v.^2, 3);
gxy = sum(abs(I_h.*I_v), 3);

% Calculate the direction of the maximum rate of change (Di Zenzo, 1986)
teta = 0.5*atan(2*gxy./(gxx - gyy));
teta_pi_half = teta+pi/2;

% F is the value of the rate of change at (x,y) for any direction
% The maximum F will be at either teta or teta+pi/2
Fteta = sqrt(0.5*((gxx+gyy) + (gxx-gyy).*cos(2*teta) + 2*gxy.*sin(2*teta)));
Fteta_pi_half = sqrt(0.5*((gxx+gyy) + (gxx-gyy).*cos(2*teta_pi_half) + 2*gxy.*sin(2*teta_pi_half)));
max_F = max([Fteta(:), Fteta_pi_half(:)],[],2);
Igrad = reshape(max_F,n_y,n_x);