function [Fn1]=ComputeF(V,x,nabla_W0,N,Xold)
Fn1=zeros(2,2,N);
nabla_W_tmp=[0,0];
Fn1_tmp=zeros(2,2);

for i=1:N
    for j=1:N
        nabla_W_tmp=nabla_W0(1:2,i,j);
        for beta=1:2
            for alpha=1:2
                uj=x(alpha,j)-Xold(alpha,j);
                ui=x(alpha,i)-Xold(alpha,i);
 Fn1_tmp(alpha,beta)= Fn1_tmp(alpha,beta)+(V(j)* (uj-ui) *nabla_W_tmp(beta));  
             end
        end
    end
    Fn1(1:2,1:2,i)=Fn1_tmp+eye(2);
    Fn1_tmp=zeros(2,2);
end
    
             

