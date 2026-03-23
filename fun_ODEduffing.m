function dydt = fun_ODEduffing(t,y,odePars,tr,fr)

    % interpolation of the applied loading
    f = interp1(tr, fr, t) ;

    % state-space equation of the Duffing oscillator
    dydt(1,1) = y(2,1) ; %dy/dt = v ;
    dydt(2,1) = odePars.m^-1 * (f - odePars.c*y(2,1) - odePars.k1*y(1,1) - odePars.k2*y(1,1)^2 - odePars.k3*y(1,1)^3) ; %dv/dt = sum(f)
end