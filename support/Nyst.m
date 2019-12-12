function [phi,Pi] = Nyst(zt,h,smp_ind)
% zt: pre-filtered image 
% h: smoothing parameter
% smp_ind: sample indices

% phi: kernel eigenvectors
% Pi: kernel eigenvalues

[M,N] = size(zt);

% parameters
ksz = 7; % Patch sizes (default 7)
krad = (ksz-1)/2;

% image padding
img = padarray(zt,[krad,krad],'symmetric');

G = fspecial('gaussian',[ksz ksz],1.2);%h2=1.2;
G = G(:);
G = G./sum(G);
Z = [];

Z = im2col(img,[ksz ksz],'sliding');
Z = Z.*repmat(G,[1 size(Z,2)]); % Gaussian weighted patches

parfor i=1:length(smp_ind)
            
        ind_y = idivide(uint32(smp_ind(i)-1), uint32(M),'floor')+1;
        ind_x = smp_ind(i)-(ind_y-1)*M;
        loc_p = [ind_x+krad,ind_y+krad];
        yp = img(loc_p(1)-krad :loc_p(1)+krad , loc_p(2)-krad :loc_p(2)+krad);
        Zc = yp(:).*G;        
        Zc = repmat(Zc,[1 M*N]);
        Ker = exp(-sum((Zc - Z).^2)/h^2);
        AB(i,:) = Ker; %[K_A,K_AB] submatices
end

K_A = AB(:,smp_ind); %K_A submatix
v = 1:M*N; 
v(smp_ind)=0; 
v= nonzeros(v);
K_AB = AB(:,v); %K_AB submatix
    
clear AB;

[phi_A,Pi] = svd(K_A); %eigenvecotrs
phi = [phi_A;K_AB'*phi_A*pinv(Pi)]; %eigenvalues

end

