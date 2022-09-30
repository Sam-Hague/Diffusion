clear; clf;

% parameters
rho = 1;     % density
cp = 1;      % specific heat
%K = 1;       % thermal conductivity

% domain and step
Lx = 10;
Ly = 10;

Nx = 51; Nt = 500;
Ny = 51;
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);              % dy should be equal to dx

% satisfy the CFL condition
c = 1;        % speed
C = 0.05;       % courant number (CFL condition C<1)
dt = C*dx/c;

% field variables
Tn = zeros(Ny,Nx);              % Temperature
x = linspace(0, Lx, Nx);       % x distance
y = linspace(0, Ly, Ny);       % y distance
[X, Y] = meshgrid(x, y);

K = ones(Ny,Nx);
K(20:25, 30:35) = 0.0001;
% Initial condition
Tn(:,:) = 0;
t = 0;

% loop
for n = 1:Nt
    tc = Tn;        % save temperature into tc for late use

    t = t+dt;       % new time
    % New temperature
    for i = 2:Nx-1
        for j = 2:Ny-1
            Tn(j,i) = tc(j,i) +...
                dt * (K(j,i)/rho/cp) *...
                ((tc(j,i+1) + tc(j+1,i) - 4*tc(j,i) + tc(j,i-1) + tc(j-1,i))/dx/dx);
        end
    end
    
    % Source term
    Sx = round(7*Nx/Nx);         % Location of source
    Sy = round(3*Ny/Ly);
    if (t<3)
        Tn(Sy,Sx) = Tn(Sy,Sx) + dt*1000/rho/cp;
    end
    % boundary conditions
    Tn(1,:) = 0;
    Tn(end,:) = 0;
    Tn(:,1) = 0;
    Tn(:,end) = Tn(:,end-1);
    
    % Visualise
    subplot(2,1,1);
    mesh(x, y, Tn); axis([0, Lx, 0, Ly, 0, 50]);
    title(sprintf('Time = %f seconds', t));
    colormap(jet)
    
    % Gradient
    [Tx, Ty] = gradient(Tn);
    qx = -K.*Tx;
    qy = -K.*Ty;
    
    subplot(2,1,2);
    imagesc(x,y,Tn);
    set(gca, 'ydir',  'norm');
    hold on
    quiver(X,Y,qx,qy,'w');
    [sx, sy] = meshgrid(2:2:8, 2:2:8);
    h = streamline(X,Y,qx,qy,sx,sy);
    set(h, 'color', 'r');
    hold off
   
    
    pause(0.001)

end
    
    
