function [energy,variance]=salient_low(img)

wavelength = 20; %10
orientation = [0 45 90 135];
g = gabor(wavelength,orientation);

gaborMag = imgaborfilt(img,g);

variance = var(gaborMag,[],3);
energy = entropy(abs(gaborMag));


