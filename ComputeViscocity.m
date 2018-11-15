function [viscosity]=ComputeViscocity(dt,eps,h,cs,V,v,alpha,Hessian_W_cor,i,j)
  viscosity=0; 
viscosity=dt*eps*h*cs*V(j)*(v(alpha,i)-v(alpha,j))*(Hessian_W_cor(1,i,j)+Hessian_W_cor(2,i,j));
    

    