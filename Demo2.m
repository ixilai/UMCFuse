close all;
clear all
clc
%%
gpuInfo = gpuDevice;
global rpos;
rpos = rotatePosition(80);
addpath function
addpath source image
%%
folder_vi = 'F:\vi_test_small\';
file_vi = dir(fullfile(folder_vi, '*.png'));

folder_ir = 'F:\ir_test_small\';
file_ir = dir(fullfile(folder_ir, '*.png'));

folder_path = 'F:\results\UMCFuse\';

vi_list= dir(strcat(folder_vi,'*.png'));
vi_num=length(vi_list);
vi_Image=cell(1,vi_num);

ir_list= dir(strcat(folder_ir,'*.png'));
ir_num=length(ir_list);
ir_Image=cell(1,ir_num);

for i=1:vi_num
    image_name = vi_list(i).name;
    image1 =  imread(strcat(folder_vi,image_name));
    image_name2 = ir_list(i).name;
    image2 =  imread(strcat(folder_ir,image_name2));
    vi_Image{i}=image1;  
    ir_Image{i}=image2;  
    fprintf('%d  %s\n',i,strcat(folder_vi,image_name));
end

j=0;
for i= 1:vi_num
    I=cell2mat(vi_Image(1,i));
    IR=cell2mat(ir_Image(1,i));

%% Add noise                                                   
sigma=0;  % noise level: 0,10,20,30,...
if sigma>0 
    v=sigma*sigma/(255*255);  
    I=imnoise(I,'gaussian',0,v);
    IR=imnoise(IR,'gaussian',0,v);
end
%% Infrared Image Pre-processing
if size(IR,3)>1
IR=rgb2gray(IR);
end
IR=im2double(IR);
tic
[~, sigma_est4, ~] = Noise_estimation(IR); %% estimate noisy level
sigma_est4 = double(sigma_est4);
y = log5(sigma_est4);
y=0.01.*round(y);
if y>0
IL = SANF(IR, y,1);
IR=IL;
end
%% Visible Image Decomposition
A1=Airlight(I,3)/255;
A=max(A1);
wsz = 3; % window size
[L1,L2,t] = tran_estimate(I);
F1=rgb2gray(I);


[~, EL1, ~] = Noise_estimation(L1); 
EL1 = double(EL1);
EL1=0.03 .*round(EL1);

k=exp(EL1);
P=0.4./exp(EL1);
P = round(P*10)/10;
%P=0.1;
E1 = SANF(L1, P ,1); 
E2 = SANF(L2, P ,1);
E3 = SANF(IR, P ,1);
D1=L1-E1;
D2=L2-E2;
D3=IR-E3;

[~, sigma_est1, ~] = Noise_estimation(D1); 
sigma_est1 = double(sigma_est1);


[~, sigma_est2, ~] = Noise_estimation(D2); 
sigma_est2 = double(sigma_est2);

[sigma_tb, sigma_est3, dis_map] = Noise_estimation(D3); 
sigma_est3 = double(sigma_est3);

e1=0.8.*log5(sigma_est1); 
e2=0.8.*log5(sigma_est2);
e3=0.8.*log5(sigma_est3);
e1=round(e1);
e2=round(e2);
e3=round(e3);

if e1>0
e11=0.1.*e1;
ED1 = SANF(D1, e11,3);
[hist_src, ~] = histcounts(D1, 256, 'Normalization', 'probability');
[hist_gamma, edges] = histcounts(ED1, 256, 'Normalization', 'probability');
kl_distance = sum((hist_src+0.1) .* log((hist_src+0.1) ./ (hist_gamma+0.1)));

if kl_distance>0.05
   ED1=D1;
end
else
    ED1=D1;
end
if e2>0
e22=0.1.*e2;
ED2 = SANF(D2, e22,3);
[hist_src2, ~] = histcounts(D2, 256, 'Normalization', 'probability');
[hist_gamma2, edges] = histcounts(ED2, 256, 'Normalization', 'probability');
kl_distance2 = sum((hist_src2+0.1) .* log((hist_src2+0.1) ./ (hist_gamma2+0.1)));

if kl_distance2>0.05 
   ED2=D2;
end
else
    ED2=D2;
end

if e3>0
e33=0.1.*e3;
ED3 = SANF(D3, e33,3);
[hist_src3, ~] = histcounts(D3, 256, 'Normalization', 'probability');
[hist_gamma3, edges] = histcounts(ED3, 256, 'Normalization', 'probability');
kl_distance3 = sum((hist_src3+0.1) .* log((hist_src3+0.1) ./ (hist_gamma3+0.1)));

if kl_distance3>0.05
   ED3=D3;
end
else
    ED3=D3;
end
%%  Low-frequency component fusion
MAP=abs(E2>E3);
E4=E2.*MAP+E3.*(1-MAP);
EN1 = entropy(E1);
EN4 = entropy(E4);
[energy1,var1]=salient_low(E1);
[energy4,var4]=salient_low(E4);

K1=exp(2.*EN1);
K4=exp(2.*EN4);
var1=mean(var1(:)).*K1;
var4=mean(var4(:)).*K4;
W3=var1./(var4+var1);
W4=var4./(var4+var1);
FE3=W4.*E4+W3.*E1;

%% High-frequency component fusion

filterbankargs.nscale = 3;
filterbankargs.norient = 1;

%   monogenic phase congruence
MPC1 = MonogenicPC(ED1, filterbankargs);
MPC2 = MonogenicPC(ED2, filterbankargs);
MPC3 = MonogenicPC(ED3, filterbankargs);
weight=MPC1+MPC2+MPC3;

weight11 = (MPC1+0.1)./max(max(weight)+0.1);
weight22 = (MPC2+0.1)./max(max(weight)+0.1);
weight33 = (MPC3+0.1)./max(max(weight)+0.1);
FD3=weight11.*ED1+weight22.*ED2+weight33.*ED3;
%% Fused Image
F=FD3+FE3;
toc
%% color
I=im2double(I);
A_YUV=ConvertRGBtoYUV(I);   
[row,column]=size(F);
F_YUV=zeros(row,column,3);
F_YUV(:,:,1)=F;
F_YUV(:,:,2)=A_YUV(:,:,2);
F_YUV(:,:,3)=A_YUV(:,:,3);
FF=ConvertYUVtoRGB(F_YUV);
%figure,imshow(FF);
%[~, name, ~] = fileparts(file_vi(i).name);
[~, name, ~] = fileparts(file_vi(i).name);
output_filename = fullfile(folder_path, [name, '.png']);
imwrite(FF, output_filename);
end

%%
function y = log5(x)
y = log(x) / log(5);
end