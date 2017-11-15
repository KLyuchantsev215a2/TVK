function [W]=ComputeW(xi,xj,h)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
% output: W = the force of the particle j effect on the initial i

W=0;
r=[0,0];
r=xi-xj;
q=norm(r,2)/h;
C=1/(3.14159265358979323846*h*h);
%if (q >=0) && (q<=3)
   % W=C*exp(-(q^2));
if (q >=0) && (q<=1)
    W=C*(10 / 7)*(1-3/2*q*q*(1-q/2));
elseif (q > 1) && (q <=2)
    W = C*(10 / 7)*(1/4)*(2 - q)^3;
end
