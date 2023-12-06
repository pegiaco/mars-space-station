function dX_dt = state_derivative_2_body(t,X)
    global mu_Mars;
    mu = mu_Mars;
    
    x = X(1,1);
    y = X(2,1);
    z = X(3,1);
    dx_dt = X(4,1);
    dy_dt = X(5,1);
    dz_dt = X(6,1);
    
    dX_dt(1:3,1) = [dx_dt, dy_dt, dz_dt];
    
    r = [x;y;z];
    dX_dt(4:6,1) = -mu/(norm(r)^3)*r;

end