function [viscocity]=ComputeViscocity(x,v,V,i,j,h,E,mu,m)
    r_ab=[0,0];
    v_ab=[0,0];
    r_ab(1)=x(1,i)-x(1,j);
    r_ab(2)=x(2,i)-x(2,j);
    v_ab(1)=v(1,i)-v(1,j);
    v_ab(2)=v(2,i)-v(2,j);
    rho_a=m/V(i);
    rho_b=m/V(j);
    
	if dot(v_ab,r_ab)>=0
		viscocity=0;
    else
        alpha=1;   
        betta=2;
        h_ab=h;
        
        ro_ab=(rho_a+rho_b)/2;
        cs_a=sqrt((E+4/3*mu)/rho_a); %скорость продольной волны
        cs_b=sqrt((E+4/3*mu)/rho_b);
        c_ab=(cs_a+ cs_b)/2;
        
        nu=h_ab*0.1;
        mu=(h_ab*dot(v_ab,r_ab)) / (dot(r_ab,r_ab) + nu^2);
        viscocity=(-alpha*c_ab*mu + betta*mu*mu) / ro_ab;
    end
    