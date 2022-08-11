% 1D Mass-Spring system with random velocity
% Written by Shin, KiBeom [irgbu1209@yonsei.ac.kr]

m = 39.948;                % Mass of atom         
N = 100;                    % Number of atom
k = 0.0509;                % Spring constat                  
a = 3.822;                 % Distance of atoms                
dt = 0.01;                 % Time                             
T = 600;                   % Tempeture       
turn = 100000;                     % Calcualated time
k_B = 8.617333262*10^(-5);          % Boltzman constant  
v = sqrt(k_B * T / m *(N-1)/N);     % Speed              

% Set velocity of each atom
velocitySeed = randn(1, N);
initVelocity = velocitySeed - mean(velocitySeed);
initVelocity = initVelocity * sqrt((N-1)*k_B*T/m) / sqrt(sum(initVelocity.^2));

% Initial condition
position = zeros(turn+1, N);
velocity = zeros(turn+1, N);
E_kin = zeros(turn+1, 1);
E_pot = zeros(turn+1, 1);

position(1, :) = a * (1:N);
distance = a * ones(1, N+1);
velocity(1, :) = initVelocity;
accel = k * diff(distance) / m;
E_kin(1, :) = m * sum(initVelocity.^2) /2;

% Main
for i = 1:turn
    tempPosition = position(i, :) + velocity(i, :) * dt + dt^2 * accel / 2;
    distance = [a diff(tempPosition) a];
    tempAccel = k * diff(distance) / m;
    tempVelocity = velocity(i, :) + dt * (tempAccel + accel) / 2;
    
    position(i+1, :) = tempPosition;
    velocity(i+1, :) = tempVelocity;
    accel = tempAccel;

    E_kin(i+1, :) = sum(tempVelocity.^2)*m/2;
    E_pot(i+1, :) = k*sum((distance-a).^2)/2;
end
E_tot = E_kin + E_pot;
t_save = 2 * E_kin / (N-1) / k_B;

% Visualization
energy = [E_kin E_pot E_tot];
t = 0:dt:turn*dt;

figure(1)  % Position & Velocity
subplot(2, 1, 1)
plot(t, position(:,1), t, position(:, 5), t, position(:,10))
xlabel('Time(ps)'), ylabel('Position(Angstrom)')
subplot(2, 1, 2)
plot(t, velocity(:, 1), t, velocity(:, 5), t, velocity(:, 10))
xlabel('Time(ps)'), ylabel('Velocity(Angstrom/ps)')

figure(2)  % Energy
plot(t, energy)
xlabel('Time(ps)'), ylabel('Energy(eV)')
legend('E_{kin}', 'E_{pot}', 'E_{tot}')

figure(3)  % Temperature
plot(t, t_save)

figure(4)  % Velocity histogram
histogram(velocity(end,:))
