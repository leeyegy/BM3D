function [V,lambda] = Orth(W_A,W_AB,m_max)
% W_A: filter submatrix
% W_AB: filter submatrix

% V: filter eigenvectors
% lambda: filter eigenvalues

W_Ah = sqrtm(pinv(W_A));
Q = W_A + W_Ah*W_AB*W_AB'*W_Ah;
[U,L,T] = svd(Q); %Eigen-decomposition of the symmetric matrix Q
V = [W_A;W_AB']*W_Ah*U*pinv(sqrt(L));% Approximated orthogonal eigenvectors
lambda = diag(L); %Approximated eigenvalues
lambda(lambda>1) = 1;

if size(V,2) > m_max 
V = V(:,1:m_max );
lambda = lambda(1:m_max);
end

end
                                                                                                                                                              