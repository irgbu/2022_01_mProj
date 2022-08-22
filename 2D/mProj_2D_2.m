% MD System with Argon Atom using Lennard-Jones Potential
% Written by Shin, KiBeom [irgbu1209@yonsei.ac.kr]

m = 39.948;                % Mass of atom         
N = 3;                     % Number of atom                 
a = 3.822;                 % Distance of atoms                
dt = 1;                    % Time                             
T = 300;                   % Tempeture
eps = 1.0325 * 10 ^(-2);   %
sig = 3.405;               %
turn = 10000;                       % Calcualated time
k_B = 8.617333262*10^(-5);          % Boltzman constant  
v = sqrt((1-1/N^2) * k_B * T / m);  % Speed 


% Set position of each atom
xPos = zeros(turn+1, N^2);
yPos = zeros(turn+1, N^2);
for i = 0:N-1
    for j = 1:N
        xPos(1, N*i+j) = a * j;
        yPos(1, N*i+j) = a * (i+1);
    end
end

xMeanPos = zeros(turn+1, 1); xMeanPos(1) = mean(xPos(1, :));
yMeanPos = zeros(turn+1, 1); yMeanPos(1) = mean(yPos(1, :));

% Set velocity of each atom
xVeloSeed = randn(1, N^2);
xInitVelo = xVeloSeed - mean(xVeloSeed);  % Set central x velocity 0

yVeloSeed = randn(1, N^2);
yInitVelo = yVeloSeed - mean(yVeloSeed);  % Set central y velocity 0

totTq = 0;
for i = 1:N^2
    dx = xPos(1, i) - xMeanPos(1);
    dy = yPos(1, i) - yMeanPos(1);
    r = sqrt(dx^2 + dy^2);
    if rem(N, 2) == 1 && i == ceil(N^2/2)
        continue
    end
    totTq = totTq + r *(yInitVelo(i)*dx/r - xInitVelo(i)*dy/r);
end

if rem(N, 2) == 0
    tq = totTq / N^2;
else
    tq = totTq / (N^2-1);
end

for i = 1:N^2
    dx = xPos(1, i) - xMeanPos(1);
    dy = yPos(1, i) - yMeanPos(1);
    r = sqrt(dx^2 + dy^2);
    if rem(N, 2) == 1 && i == ceil(N^2/2)
        continue
    end
    xInitVelo(i) = xInitVelo(i) + tq * dy / r^2;
    yInitVelo(i) = yInitVelo(i) - tq * dx / r^2;
end   % Set torque 0

xInitVelo = xInitVelo * sqrt((N^2-1)*k_B*T/m) / sqrt(sum(xInitVelo.^2));
xVelocity = zeros(turn+1, N^2); xVelocity(1, :) = xInitVelo;

yInitVelo = yInitVelo * sqrt((N^2-1)*k_B*T/m) / sqrt(sum(yInitVelo.^2));
yVelocity = zeros(turn+1, N^2); yVelocity(1, :) = yInitVelo;

xAcceleration = zeros(N);
yAcceleration = zeros(N);

E_kin = zeros(turn+1, 1);
E_kin(1, 1) = m * (sum(xInitVelo.^2)+sum(yInitVelo.^2)) / 2;
E_pot = zeros(turn+1, 1);

im = {turn};

% Main
for nTime = 1:turn
    tempE_pot = 0;
    for i = 0:N-1   % Update Position
        for j = 1:N   
            x = xPos(nTime, N*i+j);
            y = yPos(nTime, N*i+j);

            xVelo = xVelocity(nTime, N*i+j);
            yVelo = yVelocity(nTime, N*i+j);
            xAccel = xAcceleration(i+1, j);
            yAccel = yAcceleration(i+1, j);

            xPos(nTime+1, N*i+j) = x + dt*xVelo + xAccel*dt^2/2;
            yPos(nTime+1, N*i+j) = y + dt*yVelo + yAccel*dt^2/2;
        end
    end
    for i = 0:N-1   % Update Acceleration
        for j = 1:N
            % 1 - left, 2 - right, 3 - down, 4 - up
            r1 = a; r2 = a; r3 = a; r4 = a;
            a1 = 0; a2 = 0; a3 = 0; a4 = 0;
            dx1 = 0; dx2 = 0; dx3 = 0; dx4 = 0;
            dy1 = 0; dy2 = 0; dy3 = 0; dy4 = 0;
            E_pot1 = 0; E_pot2 = 0; E_pot3 = 0; E_pot4 = 0;

            xVelo = xVelocity(nTime, N*i+j);
            yVelo = yVelocity(nTime, N*i+j);
            xAccel = xAcceleration(i+1, j);
            yAccel = yAcceleration(i+1, j);
           
            x = xPos(nTime+1, N*i+j);
            y = yPos(nTime+1, N*i+j);

            if j ~= 1     % left end
                dx1 = x - xPos(nTime+1, N*i+j-1);
                dy1 = y - yPos(nTime+1, N*i+j-1);
                r1 = sqrt(dx1^2+dy1^2);
                [E_pot1,a1] = LJpotential(r1);
            end
            if j ~= N     % right end
                dx2 = x - xPos(nTime+1, N*i+j+1);
                dy2 = y - yPos(nTime+1, N*i+j+1);
                r2 = sqrt(dx2^2+dy2^2);
                [E_pot2,a2] = LJpotential(r2);
            end
            if i ~= 0     % up end
                dx3 = x - xPos(nTime+1, N*(i-1)+j);
                dy3 = y - yPos(nTime+1, N*(i-1)+j);
                r3 = sqrt(dx3^2+dy3^2);
                [E_pot3,a3] = LJpotential(r3);
            end
            if i ~= N-1   % down end
                dx4 = x - xPos(nTime+1, N*(i+1)+j);
                dy4 = y - yPos(nTime+1, N*(i+1)+j);
                r4 = sqrt(dx4^2+dy4^2);
                [E_pot4,a4] = LJpotential(r4); 
            end
            
            xAccelCache = dx1/r1*a1 + dx2/r2*a2 + dx3/r3*a3 + dx4/r4*a4;
            yAccelCache = dy1/r1*a1 + dy2/r2*a2 + dy3/r3*a3 + dy4/r4*a4;

            xVelocity(nTime+1, N*i+j) = xVelo + (xAccel+xAccelCache)*dt/2;
            yVelocity(nTime+1, N*i+j) = yVelo + (yAccel+yAccelCache)*dt/2;

            xAcceleration(i+1, j) = xAccelCache;
            yAcceleration(i+1, j) = yAccelCache;

            tempE_pot = tempE_pot + E_pot1 + E_pot2 + E_pot3 + E_pot4; 
        end
    end

    E_kin(nTime+1, :) = m * (sum(xVelocity(nTime+1, :).^2) + sum(yVelocity(nTime+1, :).^2)) / 2;
    E_pot(nTime+1, :) = tempE_pot;

    s1 = scatter(xPos(nTime, :), yPos(nTime, :), 36, [1 0 0; 0 0 1; 1 0 0; 0 0 1; 0 1 0; 0 0 1; 1 0 0; 0 0 1; 1 0 0], "filled");

    drawnow
    frame = getframe(figure(1));
    im{nTime} = frame2im(frame);

    xMeanPos(nTime+1) = mean(xPos(nTime+1, :));
    yMeanPos(nTime+1) = mean(yPos(nTime+1, :));
end

% Visualization
t = 0:dt:turn*dt;
E_tot = E_kin + E_pot;
temperature = 2 * E_kin / (N^2-1) / k_B;

figure(2)
    plot(t, E_kin, t, E_pot, t, E_tot)
    xlabel('Time(ps)'), ylabel('Energy(eV)')
    legend('E_{kin}', 'E_{pot}', 'E_{tot}')

figure(3)
    subplot(2, 1, 1)
        plot(t, xMeanPos)
        xlabel('Time(ps)'), ylabel('Position(Angstrom)')
    subplot(2, 1, 2)
        plot(t, yMeanPos)
        xlabel('Time(ps)'), ylabel('Position(Angstrom)')

figure(4)
    plot(t, temperature)



function [pot, force] = LJpotential(dist) % Lennard-Jones Potential
    global m
    eps = 1.0325 * 10 ^(-2);   %
    sig = 3.405;               %    
    pot = 4 * eps * ((sig/dist)^12 - (sig/dist)^6);
    force = (-6) * sig^6 * (dist^6 - 2 * sig^6) / dist^13 / m;
end