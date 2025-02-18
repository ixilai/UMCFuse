function t=regular(img,HazeImg,lambda)



param=0.5;
t=img;
t=im2double(t);
%%
[nRows, nCols] = size(t);
nsz = 3; NUM = nsz * nsz;

%% 算子
d{1} = [5, 5, 5; -3, 0, -3; -3, -3, -3];
d{2} = [-3, 5, 5; -3, 0, 5; -3, -3, -3];
d{3} = [-3, -3, 5; -3, 0, 5; -3, -3, 5];
d{4} = [-3, -3, -3; -3, 0, 5; -3, 5, 5];
d{5} = [5, 5, 5; -3, 0, -3; -3, -3, -3];
d{6} = [-3, -3, -3; 5, 0, -3; 5, 5, -3];
d{7} = [5, -3, -3; 5, 0, -3; 5, -3, -3];
d{8} = [5, 5, -3; 5, 0, -3; -3, -3, -3];


  %%
% normalizing filters
num_filters = length(d); 
for k = 1 : num_filters
    d{k} = d{k} / norm(d{k}(:));
end

% calculating weighting function
for k = 1 : num_filters
    WFun{k} = CalWeightFun(HazeImg, d{k}, param);
end
sobel_x=[-1 0 1; -2 0 2; -1 0 1];
sobel_y=[-1 -2 -1; 0 0 0; 1 2 1];
sobel_combined = sqrt(sobel_x.^2 + sobel_y.^2);

%%
Gx = imfilter(double(img), sobel_x);
Gy = imfilter(double(img), sobel_y);

grad_mag = sqrt(Gx.^2 + Gy.^2);
grad_mag = grad_mag / max(grad_mag(:));
weight = grad_mag.^2/ param / 2;
%% 
num_filters = length(d); 
for k = 1 : num_filters
    d{k} = d{k} / norm(d{k}(:));
end

%%
Tf = fft2(t);
DS = 0; 
for k = 1 : num_filters
    D{k} = psf2otf(d{k}, [nRows, nCols]);
    DS = DS + abs(D{k}).^2;
end
%figure,imshow(DS);
% cyclic looping for refining t
beta = 1; beta_rate = 2 * sqrt(2);
beta_max = 2^4; 
Outiter = 0;

while beta < beta_max
    % updating parameters    
    gamma = lambda / beta;
    % show the results
    Outiter = Outiter + 1; 
    %fprintf('Outer iteration %d; beta %.3g\n', Outiter, beta);    
    %figure(1000), imshow(t, []); title(num2str(Outiter)); pause(0.05);     

    % fixing t, solving u
    DU = 0;
    k=1;
    for k = 1 : num_filters
        if k<9
           dt{k} = imfilter(t, d{k}, 'circular');
           %dt{k} = imfilter(dt{k}, filter, 'replicate');
           %dt{k} = bilateralFilter2(dt{k}, dt{k}, min(dt{k}(:)), max(dt{k}(:)), 2, 0.05);
        
        else
            dt{k} = imfilter(t, d{k}, 'circular');
            %dt{k} = bilateralFilter2(dt{k}, dt{k}, min(dt{k}(:)), max(dt{k}(:)), 2, 0.5);
        end

           %dt{9} = grad_mag;
 
        u{k} = max(abs(dt{k}) - WFun{k} / beta / num_filters, 0) .* sign(dt{k}); 
        DU = DU + fft2(imfilter(u{k}, flipud(fliplr(d{k})), 'circular')); 
    end

    % fixing u, solving t;    
    t = abs(ifft2((gamma * Tf + DU) ./ ( gamma + DS)));


    % increasing beta
    beta = beta * beta_rate;
end
%%
function WFun = CalWeightFun(HazeImg, D, param)
% parameters setting
%sigma = param; 

% calculating the weighting function
HazeImg = double(HazeImg) / 255;

% weighting function
method = 'circular';
d_r = imfilter(HazeImg(:, :, 1), D, method);
d_g = imfilter(HazeImg(:, :, 2), D, method);
d_b = imfilter(HazeImg(:, :, 3), D, method);
WFun = exp(-(d_r.^2 + d_g.^2 + d_b.^2) / param / 2);
end
end
