%% Constants and Parameters
v = 0.1; %volume of parcel, m^3
l = 100;
x = 1:0.1:1000; % x (m)
t = 1:0.1:1000; % time (yrs)
dt = 0.1 %Timestep (yrs)
p = 0.3.*ones(1,numel(x));
dx = 0.2;
dtr = 0.1; %duration of rain event, minutes
I = zeros(1,numel(x));
I(1) = v.*p(1);

kfbar = 1; %Dimensionless facilitation parameter (-)
kf = kfbar/dx; %Facilitation parameter (1/m)
kcbar = 1; %Dimensionless competition parameter (-)
kc = kcbar/dx; %Competition parameter (1/m)
D  = 0.2; % Mean annual rainfall (m)
Bt = 2.1; % Transpiration rate scaling factor (-)
Bb = 0.06; % Bare soil evaporation scaling facto (-)
Bv = 0.03; % Vegetation evaporation scaling factor (-)
te = 0.97; % Annual dry duration(yr)
tr = 1-te; % Annual wet duration (yr)
P = D/tr; %Precipitation rate (m/yr)
Kmax = 0.9*P; % Maximum hydraulic conductivity (m/yr)
K0 = 0.1*P; % Hydraulic conductivity in absence of vegetation (m/yr)
Kv = Kmax - K0; %Vegetated hydraulic conductivity (m/yr)
Tmax = Bt*D/te; %Maximum transpiration rate (m/yr)
Eb = Bb*D/te; % Evaporation at bare pixel 
Ev = Bv*D/te; % Evaporation rate at vegetated pixel 
bmax = 1; %maximum biomass (-)
b = 0.2; %biomass growth per unit timestep (-)
Tcbar = 2; %Threshold transpiration rate for growth (-)
Tc = Tcbar*D/te; %Threshold transpiration rate (m/yr)

%% Create scene of vegetated (1) and bare(0) cells
b = [1 1 1 1 0 0 0 1;
    0 0 1 0 0 1 1 0;
    0 0 0 1 0 1 1 0;
    1 0 0 0 0 0 0 0;
    0 1 1 0 0 1 0 0;
    0 0 0 0 0 0 1 1;
    1 0 1 1 1 0 0 0;
    0 0 0 1 0 0 0 0];


%Create blank matrix for water balance components
Ii = zeros(size(b,1)*size(b,2),numel(t));
Ee = zeros(size(b,1)*size(b,2),numel(t));
Tt = zeros(size(b,1)*size(b,2),numel(t));

% Create blank matrix for water balance
w = zeros(size(b,1)*size(b,2),numel(t))
%% For Loop for Governing Water Balance
for j = 1:size(w,2)
    % Reshape b into vector 
    bvec =  reshape(b,1,[]);    
    %% Infitration, I
    
    for i = 2:numel(bvec)-1;
        for j = 1:numel(t)
        I(i) = I(i-1) + v.*p(i).*(1-p(i)).^(i-1)
    end
    %% Hydraulic Concuctivity, K and Transpiration,T
    
    %%%%% Ki
    % Create distance matrix
    [X,Y] = meshgrid(1:1:size(b,1));
    X = reshape(X,[],1);
    Y = reshape(Y,[],1);
    coordinates = horzcat(X,Y)
    r_mat = squareform(pdist(coordinates))
    g = repmat(reshape(b,1,[]),size(b,1)^2,1);
    g_max = ones(size(b,1)^2);
    
    % Create f function to sum
    f_real_K = g.*exp(-1*kf.^2.*(r_mat.^2));
    
    f_max_K = g_max.*exp(-1*kf.^2.*(r_mat.^2));
    f_max_vec_K = sum(f_max_K,2);
    f_max_mat_K = repmat(f_max_vec_K,1,size(b,1)^2);
    
    f_K = f_real_K./f_max_mat_K;
    sum_f_K = sum(f_K,2);
    
    % Calculate Ki from f and transform to original scene
    K = K0 + Kv.*(sum_f_K)
    %%K = vec2mat(K_vec,size(b,2))
   
    
    %%%%%% Ti
    f_real_T = g.*exp(-1*kc.^2.*(r_mat.^2));
    
    f_max_T = g_max.*exp(-1*kc.^2.*(r_mat.^2));
    f_max_vec_T = sum(f_max_T,2);
    f_max_mat_T = repmat(f_max_vec_T,1,size(b,1)^2);
    
    f_T = f_real_T./f_max_mat_T;
    sum_f_T = sum(f_T,2);
    
    % Calculate Ti from f and transform to original scene
    T = Tmax.*(sum_f_T);
    
    %%T = vec2mat(T_vec,size(b,2));
    %% Evaporation, Ei
    E = Ev.*bvec;
    indexE = find(bvec==0);
    E(indexE) = E(indexE) + Eb;
    %% Biomass Growth
    
    %b = scene
    indexb = find(T<Tc)
    b(indexb) = 0
    indexbgrowth = find(T>=Tc)
    b(indexbgrowth) = b(indexbgrowth)+0.2
    
    % Sum water balance and recreate biomass scene
    w(:,j+1) = w(:,j) + (I(:,j) - E(:,j)-T(:,j))*dt;
    b = vec2mat(b,size(b,2))
end

