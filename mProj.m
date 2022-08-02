a = 1;
k = 1;
m = 1;
T = 300;
k_B = 8.617333262*10^(-5);

% - to left, + to right

% set direction of each molcular
tempDirection = [];
for i = 1:10 
    r = randi([0 1]);
    if i > 5
        if sum(tempDirection) == 5
            r = 0;
        end
        if sum(tempDirection) + 10 - i < 5
            r = 1;
        end
    end
    direction = [tempDirection r];
    tempDirection = direction;
end

for i = 1:10
    if direction(i) == 0
        direction(i) = -1;
    end
end

v = sqrt(k_B * T);
velocity = direction * v;
coordinate = a * (1:10);
distance = a * ones([1 11]);
force = zeros([1 10]);
acceleration = force / m;
tempAcceleration = acceleration;

E_k = 10 * (1/2) * m * v^2;
E_p = 0;
E_tot = E_k + E_p;

for i = 1:100
    tempCoordinate = coordinate(i, :) + 0.1 * velocity(i, :) + (1/2) * 0.01 * acceleration(i, :); %r(t+dt)
    for j = 2:10
        distance(j) = tempCoordinate(j) - tempCoordinate(j-1);
    end % 원자간 거리 초기화
    for j = 1:10 
        tempAcceleration(j) = k * ((distance(j+1)-a)-(distance(j)-a)) / m; %a(t+dt)
    end % 가속도 초기화
    tempVelocity = velocity(i, :) + (1/2) * 0.1 * (acceleration(i, :) + tempAcceleration); %v(t+dt)
    coordinate = [coordinate;tempCoordinate]; %r(t+dt) > r(t)
    velocity = [velocity;tempVelocity]; %v(t+dt) > v(t)
    acceleration = [acceleration;tempAcceleration]; %a(t+dt) > a(t)
    
    tempE_k = (1/2) * m * sum(tempVelocity.^2);
    tempE_p = (1/2) * k * sum((distance-a).^2);
    tempE_tot = tempE_k + tempE_p;
    E_k = [E_k;tempE_k];
    E_p = [E_p;tempE_p];
    E_tot = [E_tot;tempE_tot];
end

energy = [E_k E_p E_tot];

t = 0:0.1:10;