clear; clf;

% parameters
rho = 1;
cp = 1;
K = 1;

A = K/rho/cp;

lx = 1;

nx = 21; nt = 500;
dx = lx/(nx-1);

% satisfy the CFL condition
c = 1;        % speed
C = 0.1;       % courant number (CFL condition C<1)
dt = C*dx/c;

% field variables
tn = zeros(1,nx);              % Temperature
x = linspace(0, lx, nx);       % x distance

% Initial condition
tn(:) = 0;
t = 0;
t_tot = 25;
nt = t_tot/dt;
% loop
for n = 1:nt
    tc = tn;        % save temperature into tc for late use

    t = t+dt;       % new time
    % New temperature
    for i = 2:nx-1
        tn(i) = tc(i) + dt * A * ((tc(i+1) - 2*tc(i) + tc(i-1))/dx/dx);
    end
    
    % Source term
    sx = round(7*nx/lx);         % Location of source
    %sx = 10;
    if (t<5)
        tn(sx) = tn(sx) + dt*100/rho/cp;
    end
    % boundary conditions
    %tn(1) = 0; tn(end) = 0;    % dirichlet % insulated both ends
    tn(1) = tn(2); tn(end) = 0;  % Neumann  % insulated one end
    
    % Visualise
    plot(x, tn); set(gca, 'ylim', [0, 100]);
    xlabel('Distance along rod'); ylabel('Temperature')
    title(sprintf('Time = %f seconds', t));
    pause(0.01);
end
    
    

