function [mplot]=myplot(x,V,F,N,SIG,l)

     x_coord(1:N) = x(1,1:N);
     y_coord(1:N) = x(2,1:N);
     subplot(2,2,1);
     
     scatter(x_coord,y_coord);
     axis([-l 2*l 0 2*l ]);
        
     detF=ones(1,N);
     for i = 1:N
           detF(1,i)=det(F(1:2,1:2,i));
     end
       
        tri=delaunay(x_coord,y_coord);
        subplot(2,2,2);
        trisurf(tri,x_coord,y_coord,detF(1,1:N));
        axis([0 l 0 l 0.99 1.01]);
        
        subplot(2,2,3);
        trisurf(tri,x_coord,y_coord,V(1:N,1));
        axis([0 l 0 l -0.1 0.1]);
       errSIG=ones(1,N);
       for i = 1:N
               errSIG(i)=SIG(1,1,i);
       end
       
        subplot(2,2,4);
        trisurf(tri,x_coord,y_coord,errSIG);  
        axis([0 l 0 l -1 1]);
        pause(0.0000001);
mplot=0;