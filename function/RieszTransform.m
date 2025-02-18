function [Rx, Ry] = RieszTransform(im)
% Compute the Riesz Transform of an image

im = double(im);
[rows,cols] = size(im);
[u1, u2] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
    ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)));

u1 = ifftshift(u1);
u2 = ifftshift(u2);

radius = sqrt(u1.^2 + u2.^2);
radius(1,1) = 1;

RxK = -1i*u1./radius;
RyK = -1i*u2./radius;
RxK(1,1) = 0;
RyK(1,1) = 0;

fftim = fft2(im);
Rx = real(ifft2(fftim.*RxK));
Ry = real(ifft2(fftim.*RyK));
return;
end