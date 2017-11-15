function [SIG] = ComputeStress(Fp,mu,k,N)
% Subroutine which computes Cauchy stress tensor as a function
% of the deformation gradient F
% input:  F = deformation gradient
%         mu = shear modulus
%          k = bulk modulus
% output: stress = Cauchy stresss (true stress) 
%
%   neo-Hookean material
%
SIG=zeros(2,2,N);
for i = 1:N   
         
    F = Fp(1:2,1:2,i);
    F(3,3) = 1;
    J = det(F);
    B = F*F';  % left Cauchy-Green tensor
    Biso = J^(-2/3)*B;    % isochoric part of the left Cauchy-Green tensor
    devBiso = Biso - 1/3*trace(Biso)*eye(3);   % deviatoric part Biso
    stress = (2*mu*devBiso + k/10*(J^5-J^(-5))*eye(3) )/J;
    
    SIG(1:2,1:2,i) = stress(1:2,1:2);
end