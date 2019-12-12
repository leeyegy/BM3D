function [W_A, W_AB] = Sink(phi,Pi)
% phi: kernel eigenvectors
% Pi: kernel eigenvalues

% W_A: filter submatrix
% W_AB: filter submatrix

lambda = diag(Pi);
n = size(phi,1); % number of pixels
m = size(phi,2); % number of leading eigenvectors
r = ones(n,1);

for k = 1:100
    c = 1./(phi*(lambda.*(phi'*r)));
    r = 1./(phi*(lambda.*(phi'*c)));
end

v = repmat(c,[1 m]).*phi;
parfor i = 1:m
  ABw(i,:) = r(i)*(lambda'.*phi(i,:))*v';
end
W_A = ABw(:,1:m);
W_AB = ABw(:,m+1:n);
clear ABw;
end