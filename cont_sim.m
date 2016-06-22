function dx_dot = cont_sim(t,x,f,B,zstar_init,geo_Prob,aux_Prob,Q,geo_we,geo_Ke,f_aux,df_aux,B_aux,W,dW,lambda,n,m,...
                geodesic_N,T_e,T_dot_e,geo_Aeq)
                    
z = [x(1);
       sec(x(1))*x(2);
       x(3);
       x(4)];
   
   [X, X_dot,J_opt,solved] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,zstar_init,z,T_e,T_dot_e,geo_Aeq);

    [ctrl_opt,e_rate,solved] = compute_opt_aux(aux_Prob,Q,geo_we,geo_Ke,J_opt,X,X_dot,W,dW,...
        f_aux,df_aux,B_aux,zeros(m,1),lambda,n,m);
    
    dx_dot = f(x) + B(x)*ctrl_opt;
                    
                    
end