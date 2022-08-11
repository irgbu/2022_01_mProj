% 2D Mass-Spring system
% Written by Shin, KiBeom [irgbu1209@yonsei.ac.kr]

m = 39.948;                % Mass of atom         
N = 10;                    % Number of atom
k = 0.0509;                % Spring constat                  
a = 3.822;                 % Distance of atoms                
dt = 1;                 % Time                             
T = 600;                   % Tempeture       
turn = 1000;                        % Calcualated time
k_B = 8.617333262*10^(-5);          % Boltzman constant  
v = sqrt(k_B * T / m);              % Speed     

% Set position of each atom
xPos = zeros(turn+1, 100);
yPos = zeros(turn+1, 100);
for i = 0:9
    for j = 1:10
        xPos(1, 10*i+j) = a * j;
        yPos(1, 10*i+j) = a * (i+1);
    end
end

% Set velocity of each atom
xVeloSeed = randn(1, 100);
xInitVelo = xVeloSeed - mean(xVeloSeed);
xInitVelo = xInitVelo * sqrt(100*k_B*T/m) / sqrt(sum(xInitVelo.^2));
xVelocity = zeros(turn+1, 100); xVelocity(1, :) = xInitVelo;

yVeloSeed = randn(1, 100);
yInitVelo = yVeloSeed - mean(yVeloSeed);
yInitVelo = yInitVelo * sqrt(100*k_B*T/m) / sqrt(sum(yInitVelo.^2));
yVelocity = zeros(turn+1, 100); yVelocity(1, :) = yInitVelo;

xAcceleration = zeros(10);
yAcceleration = zeros(10);

% Main
for nTime = 1:turn
    for i = 0:9         % 10
        for j = 1:10    %  1
            % 1 - left, 2 - right, 3 - down, 4 - up
            r1 = a; r2 = a; r3 = a; r4 = a;
            a1 = 0; a2 = 0; a3 = 0; a4 = 0;
            dx1 = 0; dx2 = 0; dx3 = 0; dx4 = 0;
            dy1 = 0; dy2 = 0; dy3 = 0; dy4 = 0;

            x = xPos(nTime, 10*i+j);
            y = yPos(nTime, 10*i+j);
            xVelo = xVelocity(nTime, 10*i+j);
            yVelo = yVelocity(nTime, 10*i+j);
            xAccel = xAcceleration(i+1, j);
            yAccel = yAcceleration(i+1, j);

            xPos(nTime+1, 10*i+j) = x + dt*xVelo + xAccel*dt^2/2;
            yPos(nTime+1, 10*i+j) = y + dt*yVelo + yAccel*dt^2/2;

            if j ~= 1
                dx1 = x - xPos(nTime, 10*i+j-1);
                dy1 = y - yPos(nTime, 10*i+j-1);
                r1 = sqrt(dx1^2+dy1^2);
                a1 = -k*(r1-a)/m;
            end
            if j ~= 10
                dx2 = x - xPos(nTime, 10*i+j+1);
                dy2 = y - yPos(nTime, 10*i+j+1);
                r2 = sqrt(dx2^2+dy2^2);
                a2 = -k*(r2-a)/m;
            end
            if i ~= 0
                dx3 = x - xPos(nTime, 10*(i-1)+j);
                dy3 = y - yPos(nTime, 10*(i-1)+j);
                r3 = sqrt(dx3^2+dy3^2);
                a3 = -k*(r3-a)/m;
            end
            if i ~= 9
                dx4 = x - xPos(nTime, 10*(i+1)+j);
                dy4 = y - yPos(nTime, 10*(i+1)+j);
                r4 = sqrt(dx4^2+dy4^2);
                a4 = -k*(r4-a)/m;
            end
            
            xAccelCache = dx1/r1*a1 + dx2/r2*a2 + dx3/r3*a3 + dx4/r4*a4;
            yAccelCache = dy1/r1*a1 + dy2/r2*a2 + dy3/r3*a3 + dy4/r4*a4;

            xVelocity(nTime+1, 10*i+j) = xVelo + (xAccel+xAccelCache)*dt/2;
            yVelocity(nTime+1, 10*i+j) = yVelo + (yAccel+yAccelCache)*dt/2;

            xAcceleration(i+1, j) = xAccelCache;
            yAcceleration(i+1, j) = yAccelCache;
        end
    end
    scatter(xPos(nTime, :), yPos(nTime, :))
    drawnow
    frame = getframe(figure(1));
    im{nTime} = frame2im(frame);
end
