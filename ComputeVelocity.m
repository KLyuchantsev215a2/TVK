function [velocity]=ComputeVelocity(dt,v,SIG,nabla_W,V,N,m,eps,h,Hessian_W_cor,cs)
velocity=v;
nabla_Wij=[0,0];
viscosity=0;
for i=1:N
    for j=1:N
        for beta=1:2
            for alpha=1:2
                
                nabla_Wij=nabla_W(1:2,i,j);
                viscosity=dt*eps*h*cs*V(j)*(v(alpha,i)-v(alpha,j))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));
velocity(alpha,i)=velocity(alpha,i)-dt*(V(i)*V(j)*(SIG(alpha,beta,j))*nabla_Wij(beta))/m+viscosity;
            end
        end
    end
end
