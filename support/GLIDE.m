function [zh, zt] = GLIDE(y,z,sigma)
% Global Image Denoising %
% Hossein Talebi and Peyman Milanfar, "Global Image Denoising", IEEE Transactions on Image Processing, vol 23, No. 2, pp. 755-768, February 2014.
%
% Input: y: noisy image
%        z: clean image
%        sigma: noise standard deviation
    
[M,N] = size(z);

%%%%%%%%%%% parameters %%%%%%%%%%%%%%

sd = 10; %sample distance
m0 = 5; k0 = 0; % minimum truncation and iteration
ms = 5; ks = 0.01; % truncation and iteration step sizes
m_max =200; k_max = 0.5; % minimum truncation and iteration

%%%%%%%%%%%%% Monte-Carlo %%%%%%%%%%

eps = 1; 
randn('state', 0);  % initialization
a = eps*randn(size(z(:)));
y2 = y + reshape(a,M,N);

%%%%%%%%%%%%% NLM %%%%%%%%%%%%%%%
disp(sprintf('Prefiltering by Non-Local Means...'))
tic;

%[zt,zt2] = NLM(y,y2,sigma); 
%%
[PSNR1, zt] = BM3D(z, y, sigma);
[PSNR2, zt2] = BM3D(z, y2, sigma);
%%

time = toc;
disp(sprintf('Elapsed Time = %.2f sec', time))

%%%%%%%%%%%%% GLIDE %%%%%%%%%%%%%%%
disp(sprintf('Eigen-decomposition Approximation...'))
tic;

%Sample Indices
smp_ind = Sampling(M,N,sd); %

%Adapting smoothing parameter h
h = Adapted_h(sigma);

ztt(1,:,:) = zt;
ztt(2,:,:) = zt2;

% V12 = zeros(2,M*N,m_max);
% V12 = [];
% lambda12 = zeros(m_max,1);
for k=1:2
    
zz = squeeze(ztt(k,:,:));
% Nystrom
[phi,Pi] = Nyst(zz,h,smp_ind);

% Sinkhorn
[W_A, W_AB] = Sink(phi,Pi);

% Orthogonalizing
[V12(k,:,:),lambda12(k,:)] = Orth(W_A,W_AB,m_max);

% Permutation
V12(k,:,:)  = Perm(squeeze(V12(k,:,:)),smp_ind);

end

lambda = squeeze(lambda12(1,:));
lambda2 = squeeze(lambda12(2,:));

V = squeeze(V12(1,:,:));
V2 = squeeze(V12(2,:,:));

clear V12;
clear lambda12;

time = toc;
disp(sprintf('Elapsed Time = %.2f sec', time))

% MSE minimization %
disp(sprintf('MSE minimization...'))
tic;

bd = V'*y(:);
bd2 = V2'*y2(:);

Mp=m0:ms:m_max;
Kp=k0:ks:k_max;
n = M*N;

MSE_est = zeros(length(Mp),length(Kp));

for i=1:length(Mp)
    mp = Mp(i);
    parfor j=1:length(Kp)
        
        k = Kp(j);
        lambdak = lambda(1:mp)'.^k;% diffusion
        lambdak2 = lambda2(1:mp)'.^k;% diffusion
          
        zh1 = V(:,1:mp)*(lambdak .*bd(1:mp)); 
        zh2 = V2(:,1:mp)*(lambdak2 .*bd2(1:mp));

        div = (1/(n*eps))*a'*(zh2-zh1);
        
        MSE_est(i,j) = (1/n)*sum(y(:).^2)+(1/n)* sum((lambdak.^2-2*lambdak).*bd(1:mp).^2) + 2*(sigma^2)*div - sigma^2; %SURE MSE
        
    end
end

len_m = length(m0:ms:m_max);
len_k = length(k0:ks:k_max);

[MSE_estim,I] = min(MSE_est(:));
kh = k0+ks*floor((I-1)/len_m);
mh = m0+ms*(mod(I-1,len_m));

zh = V(:,1:mh)*((lambda(1:mh)'.^kh).*(V(:,1:mh)'*y(:)));        
zh = reshape(zh,M,N);
zh(zh>255) = 255; zh(zh<0) = 0;

time = toc;
disp(sprintf('Elapsed Time = %.2f sec', time))

end
 