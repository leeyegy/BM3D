function [smp_ind] = Sampling(M,N,sd)
% M,N: Image dimensions
% sd: sample distance

% smp_ind: sample indices

c=0;
spx = [2:sd:M]; spy = [2:sd:N];
for j=1:length(spy)
    for i=1:length(spx)
        c=c+1;
        smp_ind(c) =  spx(i)+M*(spy(j)-1);
        
    end
end

smp_ind = smp_ind(smp_ind ~=1);
smp_ind = [1,smp_ind];

disp(sprintf('Sampling Percentage = %.2f %', 100*length(smp_ind)/(M*N)))

end