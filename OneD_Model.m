t = (0:0.1:100);
x = (0:0.1:100);
P = zeros(numel(t),numel(x)*100); % Plant Density
W = zeros(numel(t),numel(x)*100); % Soil Water
O = zeros(numel(t),numel(x)*100); % Surface Water

% Variable Allocation
c = 10; % conversion of uptake to growth (g/mm-m^2)
gmax = 0.05; % maximum uptake (mm/g-m^2-d)
k1 = 5; % half saturation content (mm)
Dp = 0.1; % plant dispersal (m^2/d)
a = 0.2; % max infiltration rate (1/d)
k2 = 5; % saturation constant of infiltration (g/m^2)
Wo = 0.2; % infiltration rate in absence of plant (-)
rw = 0.2; % soil water loss from ET and drainage (1/d)
Dw = 0.1; % diffusion coeff for soil water (m^2/d)
Do = 100; % diffusion coeff for surface water (m^2/d)
dt = 0.1; 
R = 5; % rainfall rate (mm/d)
d = 0.2; % loss of density from mortality (1/d)

% Initial conditions
P(1,:) = rand(1,numel(x)*100);
W(1,:) = rand(1,numel(x)*100);
O(1,:) = rand(1,numel(x)*100);

%% Hydrologically Informed Modelling from McGrath Model 

%% Update Internal Cell
for i=1:numel(t)-1;
   
    P(i+1,2:end-1) = P(i,2:end-1)...
        +(c.*gmax.*(W(i,2:end-1)./(W(i,2:end-1)+k1)).*P(i,2:end-1)-d.*P(i,2:end-1)).*dt;
        %+Dp.*((P(i,3:end)-P(i,2:end-1))./dt-(P(i,2:end-1)-P(i,2:end-2))./dt);
    W(i+1,2:end-1) = W(i,2:end-1)...
        +(a.*O(i,2:end-1).*(P(i,2:end-1)+k2.*Wo)./(P(i,2:end-1)+k2)-gmax.*(W(i,2:end-1)./(W(i,2:end-1)+k1).*P(i,2:end-1)-rw.*W(i,2:end-1)).*dt);
        %+Dw.*((W(i,3:end)-W(i,2:end-1))./dt-(W(i,2:end-1)-W(i,2:end-2))./dt);
    O(i+1,2:end-1) = O(i,2:end-1)...
        +(R - a.*O(i,2:end-1).*(P(i,2:end-1)+k2.*Wo)/(P(i,2:end-1)+k2)).*dt;
        %+Do.*((O(i,3:end)-O(i,2:end-1))./dt-(O(i,2:end-1)-O(i,2:end-2))./dt);

end
%% Separate Into Components
infiltration = a.*O.*((P+k2.*Wo)./P+k2);
infil_less_plant = a.*O;
plantuptake = gmax.*(W./(W+k1)).*P;
growth = c.*plantuptake;
dieoff = d.*P;
drainloss =rw.*W;
figure(1)
scatter(infiltration(100,:),plantuptake(100,:))
xlabel('Infiltration')
ylabel('Plant Uptake')
figure(2)
scatter(growth(100,:),dieoff(100,:))
xlabel('Growth')
ylabel('Dieoff')
figure(3)
scatter(infil_less_plant(100,:),infiltration(100,:))
xlabel('Infiltration without Plant')
ylabel('Infiltration')
figure(4)
scatter(dieoff(100,:),infiltration(100,:))
xlabel('Dieoff')
ylabel('Infiltration')
figure(5)
scatter(drainloss(100,:),plantuptake(100,:))
xlabel('Drainage loss')
ylabel('Plant Uptake')
%% convert to x and y 

Pxy100 = vec2mat(P(100,:),100);
Pxy200 = vec2mat(P(200,:),100);
Pxy300 = vec2mat(P(300,:),100);
Pxy400 = vec2mat(P(400,:),100);
Pxy800 = vec2mat(P(800,:),100);
mesh(Pxy100)
mesh(Pxy200)
mesh(Pxy300)
mesh(Pxy800)


%% Boundary Conditions
scatter(P(100,:),O(100,:));
z = mean(P,2);
scatter(mean(W,2),mean(P,2));
%%scatter(P(:,1),W(:,1))
%%mesh(P)
