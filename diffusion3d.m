clear; clf;

% parameters
rho = 1.205;     % density
cp = 1005;      % specific heat
%K = 1;       % thermal conductivity

% domain and step
Lx = 10;
Ly = 10;
Lz = 10;

Nx = 21; Nt = 120;
Ny = 21; 
Nz = 21;
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);              % dy should be equal to dx
dz = Lz/(Nz-1);

% satisfy the CFL condition
c = 1;        % speed
C = 0.05;       % courant number (CFL condition C<1)
dt = C*dx/c;
dt = 10;

% field variables
Tn = zeros(Nx, Ny, Nz);        % Temperature
x = linspace(0, Lx, Nx);       % x distance
y = linspace(0, Ly, Ny);       % y distance
z = linspace(0, Lz, Nz);       % z distance
[X, Y, Z] = meshgrid(x, y, z);

K = ones(Nx,Ny,Nz) + 0.0257*100;
K([1 end],:,:) = 0.001;
K(:,[1 end], :) = 0.001;
K(:,:,[1 end]) = 0.001;
% Initial condition
Tn(:,:,:) = 25;
t = 0;

% loop
for n = 1:Nt
    Tc = Tn;        % save temperature into tc for late use

    t = t+dt;       % new time
    % New temperature
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 2:Nz-1
                Tn(i, j, k) = Tc(i, j, k) +...
                    dt * (K(i,j,k)/rho/cp) *...
                    ((Tc(i+1,j,k) - 2*Tc(i,j,k) + Tc(i-1,j,k))/dx/dx) + ...
                    ((Tc(i,j+1,k) - 2*Tc(i,j,k) + Tc(i,j-1,k))/dy/dy) + ...
                    ((Tc(i,j,k+1) - 2*Tc(i,j,k) + Tc(i,j,k-1))/dz/dz);
            end
        end
    end
    
    Tbar = mean(Tn(:));
    
    % Source term
    if (t < 10*3600)
        Tn(10,10,2) = Tn(10,10,2) + dt*1000/rho/cp;
    end
    % boundary conditions
    Tn(1,:,:) = Tn(2,:,:);
    Tn(end,:,:) = Tn(end-1,:,:);
    Tn(:,1,:) = Tn(:,2,:);
    Tn(:,end,:) = Tn(:,end-1,:);
    Tn(:,:,1) = Tn(:,:,2);
    Tn(:,:,end) = Tn(:,:,end-1);
    
    % Sink term
    if (t > 5*3600)
        %Tn(10,10,2) = Tn(10,10,2) + dt*1000/rho/cp;
        Tn(end,9:11,9:11) = 20;
    end
    
    if (mod(t, 600) == 0)
        clf;
        % Visualise
        subplot(2,1,1);
        slice(X,Y,Z,Tn,5,5,2); colorbar;
        axis([0 Lx 0 Ly 0 Lz]);
        title(sprintf('Average temperature = %.2f Time = %f minutes', Tbar, t/60));

    
        % Gradient
        [Tx, Ty, Tz] = gradient(Tn);
        qx = -K.*Tx;
        qy = -K.*Ty;
        qz = -K.*Tz;

        subplot(2,1,2);
        streamslice(X,Y,Z,qx,qy,qz,9.9,9.9,0.1);
        axis([0 Lx 0 Ly 0 Lz]);
        hold on

        imagesc(x,y,Tn);
        set(gca, 'ydir',  'norm');
        hold on

        [sx, sy] = meshgrid([2,5,8], [2,5,8], [2,5,8]);
        h = streamline(X,Y,Z,qx,qy,qz,sx,sy,sz);
        set(h, 'color', 'r');
        hold off
   
    pause(0.001)
    end
end
    
    
