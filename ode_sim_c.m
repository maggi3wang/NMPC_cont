function dx_dot = ode_sim_c(t,x,f,B,...
  geo_Prob,W,m,n,geo_N,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq,rho)

[X, X_dot,J_opt,exitflag] = compute_geodesic_tom...
  (geo_Prob,n,geo_N,zeros(n,1),x,T_e,T_dot_e,geo_Aeq);

u = aux_control(rho,X,X_dot,geo_Ke,geo_we,...
        B,W,m);

dx_dot = f(x) + B*u;

end