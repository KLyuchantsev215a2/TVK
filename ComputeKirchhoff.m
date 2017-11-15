function [P] = ComputeKirchhoff(F,SIG,N)
P=zeros(2,2,N);
for i = 1:N   
    P(1:2,1:2,i)=F(1:2,1:2,i)*SIG(1:2,1:2,i);
end