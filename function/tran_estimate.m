function [L1,L2,t4]=tran_estimate(I)
A1=Airlight(I,3)/255;
A=max(A1);

I=im2double(I);

min_I=min(I,[],3);
%%
B=1.2778;

min_J4=1./(log(A-min_I).^B);
t4=exp(log(A-min_I)./(1-min_J4./A));
t4=real(t4);
%%
se = strel('square', 3);
t4=imclose(t4, se); 
%%
[~, lambda, ~] = Noise_estimation(t4);
lambda = double(lambda);
lambda=1.5.*exp(-lambda);
%%
t4=regular(t4,I,lambda);
f1=rgb2gray(I);
L1=t4.*f1;
L2=(1-t4).*f1;

