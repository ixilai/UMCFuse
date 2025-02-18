%% The implementation of Scale-Aware Edge-Preserving (SAEG) filter
% I — Input image; 
% lambda — Global smoothness weight;
% r1 — Radius of local window (i.e., the scale parameter),
%       The information whose scale is smaller than r1 would be smoothed
% r0 — The mimimum scale considerded in the filtering,
%       i.e., the scales within the range [r0, r1] are smoothed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = SANF(I, lambda, r1, r0)
    if (~exist('r0','var'))
       r0 = 1;
    end
    S = I;
    ch = size(I,3);
    for c = 1:ch
        wtbi = ones(size(I(:,:,c)));
        for sigma = r0:r1
            wtbo = OriGrad_MD(I(:,:,c), sigma);
            wtbi = min(wtbo, wtbi);
        end
        wtbi = medfilt2(wtbi, [r1,r1]);  % Applying median filter(optionally)
        S(:,:,c) = salient(I(:,:,c), lambda, wtbi);
    end

end

function wtb = OriGrad_MD(fin, sigma)
    global rpos;
    vareps_s =20e-3; %1e-5
    fin0 = fin;
%    fin = log(fin+eps);
    sigma2 = min(sigma, 50);
%     if sigma2<4
%         sigma2 = 2;
%     else
       sigma2 = 4*sigma*10; 
%     end
   if sigma2~=0
       if sigma<4
           if sigma<=1
               w = floor(1*sigma); %w = floor(1*sigma-1)
           else
                w = floor(1*sigma); %w = floor(0.5*sigma-1)
           end
       else
           w = floor(1*sigma);
       end

       h = fspecial('gaussian', [2*w+1, 2*w+1], sigma2);
       sfin = imfilter(fin, h, 'symmetric');
   else
       sfin = fin;
   end
   dx = diff(sfin,1,2);
   dx = padarray(dx, [0 1 0], 'post');
   dy = diff(sfin,1,1);
   dy = padarray(dy, [1 0 0], 'post');
%     [dx,dy] = GradientMethod(sfin, 'sobel');
   
%    sigma1 = 0;
%    w = floor(max(1*sigma1, 0));
%    h = fspecial('gaussian', [2*w+1, 2*w+1], sigma1+1);
%    dxx = imfilter(dx.*dx, h);
%    dxy = imfilter(dx.*dy, h);
%    dyy = imfilter(dy.*dy, h);   
% 
%     x1=dxy;
%     leda = sqrt( (dxx-dyy).^2+4*dxy.^2 );
%     leda = (dxx+dyy+leda)*0.5;
%     x2 = leda - dxx;
%     dxy0 = (dxy==0);
%     ps = dxx(dxy0) > dyy(dxy0);
%     tt = x1(dxy0); tt(ps) = 1; x1(dxy0) = tt;
%     tt = x2(dxy0); tt(ps) = 0; x2(dxy0) = tt;    
%     ps = dxx(dxy0) < dyy(dxy0);
%     tt = x1(dxy0); tt(ps) = 0; x1(dxy0) = tt;
%     tt = x2(dxy0); tt(ps) = 1; x2(dxy0) = tt;   
%     ps = dxx(dxy0) == dyy(dxy0);
%     tt = x1(dxy0); tt(ps) = 1; x1(dxy0) = tt;
%     tt = x2(dxy0); tt(ps) = 1; x2(dxy0) = tt;
    
    x1 = dx; x2 = dy;
    mag = sqrt(x1.^2 + x2.^2);
    x1 = x1./mag;
    x2 = x2./mag;
    
%     x1 = ones(size(x1))*cosd(90); x2 = ones(size(x2))*sind(90);
    
%     para.fig = 'Example';
%     mainori = x1 + x2*1i;
%     ShowImageGrad(fin0, para, 20*mainori);

    dg = atand(x2./x1); 
    npt = dg<0; dg(npt) = dg(npt)+180;
    ndg = floor(2*dg+0.5);

   sigma0 = sigma*10;
%    if sigma0<4
%        sigma0 = 0.25;
%    else
       sigma0 = max(sigma, 0.25);
%    end
   if sigma0~=0
       w = floor(0.25*sigma); 
%        w = ceil(0.5*sigma); 
       h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);
       ssfin = imfilter(fin, h, 'symmetric');
   else
       ssfin = fin;
   end
   dx = diff(ssfin,1,2);
   dx = padarray(dx, [0 1 0], 'post');
   dy = diff(ssfin,1,1);
   dy = padarray(dy, [1 0 0], 'post');
     
   wi = ceil(sigma);
%      [AS, SA] = analyGrad_1DG(dx, dy, x1, x2, wi);
%    [AS, SA] = analyGrad_1DG_MD(ssfin, ndg, wi);
%    [AS, SA] = analyGrad_1DG_MD1(dx, dy, x1, x2, ndg, wi);
   [AS, SA] = analyGrad_C(ssfin, ndg, rpos, wi);
%    [AS, SA] = analyGrad_CINT(ssfin, ndg, rpos, wi);
    
    r=6;
    AS(:, 1:r) = AS(:, r+1:2*r);  SA(:, 1:r) = SA(:, r+1:2*r);  
   
   wtb = AS./max(SA, vareps_s);
%    wtb = exp(-(1-wtb).^2/(2*0.05*0.05));
end