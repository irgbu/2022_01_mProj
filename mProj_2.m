% 1D Mass-Spring system with random velocity
% Written by Shin, KiBeom [irgbu1209@yonsei.ac.kr]

m = 1;                     % Mass of atom         
N = 10000;                 % Number of atom
k = 1;                     % Spring constat                  
a = 1;                     % Distance of atoms                
dt = 0.01;                 % Time                             
T = 600;                   % Tempeture       
turn = 10000;                       % Calcualated time
k_B = 8.617333262*10^(-5);          % Boltzman constant  
v = sqrt(k_B * T / m *(N-1)/N);     % Speed

% Set direction of each atom
direction = ones(1, N);
directionCache = randperm(N, N/2);
for i = 1:N/2
    direction(directionCache(i)) = -1;
end

% Initial condition
coordi = zeros(turn+1, N);
velo = zeros(turn+1, N);
E_kin = zeros(turn+1, 1);
E_pot = zeros(turn+1, 1);

coordi(1, :) = a * (1:N);
distance = a * ones(1, N+1);
velo(1, :) = v * direction;
accel = k * diff(distance) / m;
E_kin(1, :) = N * m * v^2 /2;

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
end
E_tot = E_kin + E_pot;
t_save = 2 * E_kin / (N-1) / k_B;

% Visualization
energy = [E_kin E_pot E_tot];
t = 0:dt:turn*dt;

figure(1)  % Position & Velocity
subplot(2, 1, 1)
plot(t, coordi(:, 1), t, coordi(:, 5), t, coordi(:, 10))
xlabel('Time(ps)'), ylabel('Position(Angstrom)')
subplot(2, 1, 2)
plot(t, velo(:, 1), t, velo(:, 5), t, velo(:, 10))
xlabel('Time(ps)'), ylabel('Velocity(Angstrom/ps)')

figure(2)  % Energy
plot(t, energy)
xlabel('Time(ps)'), ylabel('Energy(eV)')
legend('E_{kin}', 'E_{pot}', 'E_{tot}')

figure(3)  % Temperature
plot(t, t_save)

figure(4)  % Velocity histogram
histogram(velo(end,:))
