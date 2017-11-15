 function [V] = computeV(N,W_cor,m)

V=zeros(N,1);
rho=zeros(N,1);
for i = 1:N
    for j = 1:N
       rho(i)=rho(i)+m*W_cor(i,j); 
    end
end

for i = 1:N
    V(i)=m/rho(i);
end
