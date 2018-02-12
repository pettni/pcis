function con = constants
    % Common parameters
    con.m = 1800; %kg

    % ACC parameters
    con.f0bar = 74.63;
    con.f1bar = 40.59;

    % ACC bounds
    con.u_min = 25/3.6;
    con.u_max = 30/3.6;
    con.h_min = 5;
    con.vl = 7;

    con.Fw_max = 2*con.m;
    con.Fw_min = -3*con.m;

    % LK parameters
    con.Iz = 3270; % kgm^2
    con.a = 1.2; % m
    con.b = 1.65; % m
    con.Caf = 1.4e5; % N/rad 
    con.Car = 1.2e5; % N/rad

    % LK bounds
    con.y_max = 0.7;
    con.nu_max = 1;             % max lateral speed [m/s]
    con.psi_max = deg2rad(20);  % max yaw angle     [deg]
    con.r_max = deg2rad(40);    % max yaw rate      [deg/s]
    
    con.df_max = deg2rad(20);        % control
    con.rd_max = 0.06 * con.u_min;   % inverse curvature at lowest speed
    con.rd_min = 0.06 * con.u_max;   % inverse curvature at highest speed

    % Time discretization
    con.dt = 0.1;

    % Rho factors for convergence
    con.rho_lk = 0.02*[con.y_max; con.nu_max; con.psi_max; con.r_max];  % finishes at iteration 16 (39000 seconds)
    con.rho_acc = [0.05; 0.05];
end
