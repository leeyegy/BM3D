 function [Vnew] = Perm(V,indx)
 % V: filter eigenvectors
% indx: sample indices

% V: permuted filter eigenvectors


     c=0;
     Vnew = V;
     
     if indx(end) == size(V,1)
         q = length(indx)-1;
         indxx = indx;
     else
         q = length(indx);
         indxx = [indx,size(V,1)];
     end
     
     in = length(indx);

for i=1:q
            c=c+1;

            Vnew(indxx(c),:) = V(c,:);            
            l = indxx(c+1)-indxx(c)-1;
            t1 = indxx(c)+1; t2 = indxx(c+1)-1; dt(i) = t2-t1;
            Vnew(t1:t2,:) = V(in+1:in+l,:);
            in = in + l;
                          
end

if indx(end) == size(V,1);
         Vnew(end,:) = V(c+1,:); 
end

clear V;

 end
 
 
