 function [zt1,zt2] = NLM(y1,y2,sigma)

% y1: input noisy image 
% y2: input noisy image 
% sigma: standard deviation of noise

% zt1: denoised y1
% zt2: denoised y2

[M,N] = size(y1);

%%%%%%%%%%%%%% smoothing parameter & patch size %%%%%%%%%%%%%%%
% psz: patch size
% h:smoothing parameter
if sigma<5 
    h=0.9*sigma; psz=3;
elseif sigma<15
    h=0.85*sigma; psz=3;
elseif sigma<25
    h=0.75*sigma; psz=5;
elseif sigma<35
    h=0.65*sigma; psz=5;
else h=0.55*sigma; psz=5;
end

R=10; % Search radius
wsz = 2*R+1; % Search window
prad = (psz-1)/2;

%%%%%%%%%%%%% padding %%%%%%%%%%%%%%
Img = padarray(y1,[R+prad R+prad],'symmetric');
Img(:,:,2) = padarray(y2,[R+prad R+prad],'symmetric');

Z = [];

for k=1:2
    img = squeeze(Img(:,:,k));
  for i=1:M
    
    parfor j=1:N
                
    loc_p = [i+R+prad,j+R+prad];
    
    yp = img(loc_p(1)-R :loc_p(1)+R , loc_p(2)-R :loc_p(2)+R);

    
    rowMin = loc_p(1,1) - R - prad;
    rowMax = loc_p(1,1) + R + prad;
    colMin = loc_p(1,2) - R - prad;
    colMax = loc_p(1,2) + R + prad;
    
    z = img(rowMin:rowMax,colMin:colMax);
    Z = im2col(z,[psz psz],'sliding');
    
    Zc = Z(:,(wsz^2-1)/2+1);
    Zc = repmat(Zc,[1 wsz^2]);    
    d2 =  mean((Zc - Z).^2);
    Ker = exp(-max(d2-2*sigma^2,0)/h^2);

    Ker = Ker/sum(Ker);
    zt(i,j,k) =Ker*yp(:);
    
    end
    
  end
  
end

zt(zt>255) = 255; zt(zt<0) = 0;

zt1 = squeeze(zt(:,:,1));
zt2 = squeeze(zt(:,:,2));

 end

