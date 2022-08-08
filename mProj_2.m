% 1D Mass-Spring system with real values
% Made by Shin, KiBeom [irgbu1209@yonsei.ac.kr]

m = 39.948;         % Mass of atom         
N = 10;        % Number of atom
k = 0.0509;         % Spring constat                  
a = 3.822;         % Distance of atoms                
dt = 0.1;      % Time                             
T = 300;       % Tempeture       
turn = 1000;   % Calcualated time
k_B = 8.617333262*10^(-5); % Boltzman constant  
v = sqrt(k_B * T / m);     % Speed              

% Set direction of each atom
direction = zeros(1, 10);
for i = 1:10
    r = randi([0 1]);
    if i > 5
        if sum(direction) == 5
            r = 0;
        end
        if sum(direction) + 10 - i < 5
            r = 1;
        end
    end
    direction(i) = r;
end
for i = 1:10
    if direction(i) == 0
        direction(i) = -1;
    end
end

% Initial condition
coordi = zeros(turn+1, N);
velo = zeros(turn+1, N);
E_kin = zeros(turn+1, 1);
E_pot = zeros(turn+1, 1);
E_tot = zeros(turn+1, 1);

coordi(1, :) = a * (1:N);                
distance = a * ones(1, 11);              
velo(1, :) = v * direction;              
accel = k * diff(distance) / m;          
E_kin(1, :) = 10 * m * v^2 /2;         
E_tot(1, :) = E_kin(1, :) + E_pot(1, :);

for i = 1:turn
    tempCoordi = coordi(i, :) + velo(i, :) * dt + dt^2 * accel / 2;
    distance = [a diff(tempCoordi) a];
    tempAccel = k * diff(distance) / m;
    tempVelo = velo(i, :) + dt * (tempAccel + accel) / 2;
    
    coordi(i+1, :) = tempCoordi;
    velo(i+1, :) = tempVelo;
    accel = tempAccel;

    E_kin(i+1, :) = sum(tempVelo.^2)*m/2;
    E_pot(i+1, :) = k*sum((distance-a).^2)/2;
    E_tot(i+1, :) = E_kin(i+1)+E_pot(i+1);
end

% Visualization
energy = [E_kin E_pot E_tot];
t = 0:dt:turn*dt;

figure(1)  % Position & Velocity
subplot(2, 1, 1)
plot(t, coordi(:, 1), t, coordi(:, 5), t, coordi(:, 10))
xlabel('Time'), ylabel('Position')
subplot(2, 1, 2)
plot(t, velo(:, 1), t, velo(:, 5), t, velo(:, 10))
xlabel('Time'), ylabel('Velocity')

figure(2)  % Energy
plot(t, energy)
xlabel('Time'), ylabel('Energy')
legend('E_{kin}', 'E_{pot}', 'E_{tot}')
