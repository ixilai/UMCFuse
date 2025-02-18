function [Gxy]=gradientxy(image)
% Caculates the gradient magnitude of an image.

horizontalSobel=[-3,0,3; -10,0,10; -3,0,3];
verticalSobel=horizontalSobel';
Gx = filter2(horizontalSobel,image);
Gy = filter2(verticalSobel,image);
Gxy=sqrt(Gx.*Gx+Gy.*Gy);
return;
end