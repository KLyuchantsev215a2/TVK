function [Ln1]=ComputeL(v,V,nabla_W_cor,N)

L=zeros(2,2,N);
nabla_W_tmp=[0,0];
L_tmp=zeros(2,2);
v_tmp=[0,0];


for i = 1:N
    for j = 1:N        
         v_tmpj=v(1:2,j);
         
         nabla_W_tmp=nabla_W_cor(1:2,i,j); 
        
         for beta=1:2   
             for alpha=1:2
                 L_tmp(alpha,beta)=L_tmp(alpha,beta)+(V(j))*(v_tmpj(alpha))*nabla_W_tmp(beta); 
             end
         end
        
    end
    L(1:2,1:2,i)=L_tmp;
    L_tmp=zeros(2,2);
end

Ln1=L;