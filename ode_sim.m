function dx_dot = ode_sim(t,x,u,f,B,B_w,w_max)

dx_dot = f(x) + B*u + B_w*w_max;%*cos((2*pi/1)*t + 0.13);

end