function [W_cor]=ComputeW_cor(N,x,xper,V,h)

% for i = 1:N  
%     xi=xper(1:2,i);
%            for j = 1:N
%                xj=x(1:2,j);
%                 W_cor(i,j)=ComputeW(xi,xj,h);
%            end
% end
% return;

xi=[0,0];
xj=[0,0];

for i = 1:N
    sumW=zeros(2,1);
    alpha=0;
    betaij=zeros(2,2);
    cormat=zeros(2,2);
    xi=xper(1:2,i);
           for j = 1:N
                   xj=x(1:2,j);
                   ri=xi-xj;
                   sumW=sumW+V(j)*ComputeW(xi,xj,h)*ri;
            end
            
            for beta=1:2
                for j = 1:N
                    xj=x(1:2,j);
                    ri=xi-xj;
                    cormat=cormat+(ri*ri'*V(j)*ComputeW(xi,xj,h));%51
                 end
            end
            
            betaij=cormat^(-1)*sumW;
            
            for j = 1:N
                 xj=x(1:2,j);
                 ri=xi-xj;
                 alpha=alpha+(V(j)*(1+dot(betaij,ri))*ComputeW(xi,xj,h));%52
            end
            alpha=1/alpha;
            
            for j = 1:N
                 xj=x(1:2,j);
                 ri=xi-xj;
                 W_cor(i,j)=ComputeW(xi,xj,h)*alpha*(1+dot(betaij,ri));%47
            end
       
end


