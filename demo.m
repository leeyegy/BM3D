% GLIDE MATLAB toolbox

% Hossein Talebi and Peyman Milanfar, "Global Image Denoising", IEEE Transactions on Image Processing, vol 23, No. 2, pp. 755-768, February 2014.
% This is experimental software. It is provided for non-commercial research purposes only. Use at your own risk. No warranty is implied by 
% this distribution. Copyright ?2014 by University of California.
addpath('BM3D_images');
addpath('support');
clc;  
%z = double(imread('monarch.png')); % clean image
z = double(imread('Cameraman256.png')); % clean image
sigma = 50;
randn('state', 1); % initialization
y = z + randn(size(z)) * sigma; % noisy image

% matlabpool(2) % Uncomment for parallel computation
[zh, zt] = GLIDE(y,z,sigma);

% matlabpool close

% Display Noisy Image
en = y - z;
MSE_Noisy = mean(en(:).^2);
PSNR_Noisy = 10*log10((255^2)/MSE_Noisy);
figure,imshow(y,[0 255]);
title(sprintf('Noisy Image (PSNR = %.2f dB)', PSNR_Noisy), 'FontSize', 12);

% Display Prefiltered Image
et = zt - z;
MSE_PreFilter = mean(et(:).^2);
PSNR_PreFilter = 10*log10((255^2)/MSE_PreFilter);
figure,imshow(zt,[0 255]);
title(sprintf('PreFiltered Image (PSNR = %.2f dB)', PSNR_PreFilter), 'FontSize', 12);

% Display GLIDE Image
eh = zh - z;
MSE_GLIDE = mean(eh(:).^2);
PSNR_GLIDE = 10*log10((255^2)/MSE_GLIDE);
disp(sprintf('GLIDE PSNR = %.2f dB', PSNR_GLIDE))
figure,imshow(zh,[0 255]);
title(sprintf('GLIDE (PSNR = %.2f dB)', PSNR_GLIDE), 'FontSize', 12);
