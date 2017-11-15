function [f] = ComputeForse(V,P,nabla_W_cor,N)
f=zeros(2,N);
for i = 1:N  
    for j=1:N
        for beta=1:2
            for alpha=1:2
             f(alpha,i)=f(alpha,i)+V(i)*V(j)*(P(alpha,beta,i)*nabla_W_cor(beta,i)-P(alpha,beta,j)*nabla_W_cor(beta,j));
            end
        end
    end
end