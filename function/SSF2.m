
function S=SSF2(I,I0)
global rpos;
rpos = rotatePosition(80);

[N,M,D]  = size(I);
sizeI2D  = [N,M];  

Dx       = [1, -1]/2;
Dy       = Dx';
otfDx    = psf2otf(Dx, sizeI2D); 
otfDy    = psf2otf(Dy, sizeI2D);

fxx      = [1, -2, 1]/4;
fyy      = fxx';
fuu      = [1 0  0;  0 -2 0;  0 0 1]/4;
fvv      = [0 0  1;  0 -2 0;  1 0 0]/4;
otfFxx   = psf2otf(fxx, sizeI2D);
otfFyy   = psf2otf(fyy, sizeI2D);
otfFuu   = psf2otf(fuu, sizeI2D);
otfFvv   = psf2otf(fvv, sizeI2D);

Denormin1 = abs(otfDx).^2 + abs(otfDy).^2;
Denormin2 = abs(otfFxx).^2 + abs(otfFyy).^2+abs(otfFuu).^2 + abs(otfFvv).^2;
if D>1
    Denormin1 = repmat(Denormin1,[1,1,D]);
    Denormin2 = repmat(Denormin2,[1,1,D]);
end

S           = I;
iter        = 1;
errs        = []; 
alpha       = 0.1;
beta        = 0.02;
kappa       = 1.2;
tau         = 0.95;
iter_max    = 40;
lambda_max  = 1e8;
lambda      = 10*beta;
Normin0     = fft2(S);
fprintf('The code is runing and may take a few minutes ... \n');
while (lambda <= lambda_max) && (iter <= iter_max)

    Denormin  = 1.0 + alpha*Denormin1 + lambda*Denormin2;
    
    % first-order gradients
    gx = imfilter(S, Dx, 'circular');
    gy = imfilter(S, Dy, 'circular');
    
    % second-order gradients
    gxx = imfilter(S, fxx, 'circular');
    gyy = imfilter(S, fyy, 'circular');
    guu = imfilter(S, fuu, 'circular');
    gvv = imfilter(S, fvv, 'circular');
   
    if D==1
        t = (gxx.^2+gyy.^2+guu.^2+gvv.^2)<beta/lambda;
    else
        t = sum((gxx.^2+gyy.^2+guu.^2+gvv.^2),3)<beta/lambda;
        t = repmat(t,[1,1,D]);
    end
    gxx(t)=0; gyy(t)=0; guu(t)=0; gvv(t)=0;
    
    
    Normin_x = circshift(imfilter(gx, Dx(end:-1:1), 'circular'),[0, 1]);
    Normin_y = circshift(imfilter(gy, Dy(end:-1:1), 'circular'),[1, 0]);
    Normin1 = Normin_x+Normin_y;
    
    Normin_xx = imfilter(gxx, fxx(end:-1:1,end:-1:1), 'circular');
    Normin_yy = imfilter(gyy, fyy(end:-1:1,end:-1:1), 'circular');
    Normin_uu = imfilter(guu, fuu(end:-1:1,end:-1:1), 'circular');
    Normin_vv = imfilter(gvv, fvv(end:-1:1,end:-1:1), 'circular');
    
    Normin2 = Normin_xx+Normin_yy+Normin_uu+Normin_vv;
    FS = (Normin0 + alpha*fft2(Normin1) + lambda*fft2(Normin2))./Denormin;
    S  = real(ifft2(FS));
    
    err_vec = (S-I0).^2;
    errs(iter)= mse(err_vec(:));

    alpha  = tau*alpha;
    lambda = kappa*lambda;
    iter = iter+1;
    if (mod(iter,50))
        fprintf('.');
    else
        fprintf('\n');
    end
end
fprintf('\n');