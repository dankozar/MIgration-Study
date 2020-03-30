%% Constants and Parameters
clear, clc
v = 0.1; %volume of parcel, m^2 
nt = 1:0.1:1000; % time (yrs)
dt = 0.1; %Timestep (yrs)
dx = 5; % grid size (m)

kfbar = 1; %Dimensionless facilitation parameter (-)
kf = kfbar/dx; %Facilitation parameter (1/m)
kcbar = 1; %Dimensionless competition parameter (-)
kc = kcbar/dx; %Competition parameter (1/m)
D  = 0.2; % Mean annual rainfall (m)
Bt = 2.1; % Transpiration rate scaling factor (-)
Bb = 0.06; % Bare soil evaporation scaling factor (-)
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
Tcbar = 2; %Threshold transpiration rate for growth (-)
Tc = Tcbar*D/te; %Threshold transpiration rate (m/yr)
wcbar = 0.048; % Dimensionless water threshold for plant growth (-)
wc = wcbar*dx; % Threshold water for biomass growth in bare area (m)
growth = 0.2; % Biomass growth increment (kg/m2-yr)
blimit = 2; % Biomass density limit (kg/m2-yr)
%% Create scene of vegetated (1) and bare(0) cells
b = [1 0 0;
     0 0 1;
     0 0 0];
% uniform probablity at t=0 
p = 0.3.*ones(size(b,1),size(b,2));

%Create blank matrix for water balance components
Imat = zeros(size(b,1)*size(b,2),numel(nt));
Emat = zeros(size(b,1)*size(b,2),numel(nt));
Tmat = zeros(size(b,1)*size(b,2),numel(nt));
Kmat = zeros(size(b,1)*size(b,2),numel(nt));
bmat = zeros(size(b,1)*size(b,2),numel(nt));
bvecmat = zeros(size(b,1)*size(b,2),numel(nt));

% Create blank matrix for water balance
w = zeros(size(b,1)*size(b,2),numel(nt));
x = randn(1,numel(b)); 
w(:,1) = abs(x)./max(abs(x));
%% For Loop for Governing Water Balance
for t=1:numel(nt);
   %% create b vector 
   bvec = zeros(size(b,1)*size(b,2),1); % blank space for b vector
   bvec(b ~= 0) = 1; % if the value of b is not zero, assign the cell a 1 in bvec
   bvecmat(:,t) = bvec; % save bvec in matrix for each time step
    %% Infitration, I
    I = zeros(size(b,1),size(b,2)); % empty I space
    runoff = zeros(size(b,1),size(b,2)); % empty runoff space
    for i = 2:size(b,1);
        I(i,:) = I(i-1,:) + v.*p(i,:).*(1-p(i,:)).^(i-1);
        runoff(i,:) = v+(v-I(i-1,:))-I(i,:);
    end
    
    I(1,:) = (v+runoff(end,:)).*(1-p(1,:)); % assign first cell in each column with value from last row
    I = reshape(I,1,[]); % reshape I into vector 
    
    Imat(:,t) = I; % save I for each time step
    %% Hydraulic Concuctivity, K and Transpiration,T
    
    %%%%% Ki
    % Create distance matrix
    [X,Y] = meshgrid(1:1:size(b,1)); % create meshgrid of same dimensions as b
    X = reshape(X,[],1); % save x values for each cell
    Y = reshape(Y,[],1); % savey values for each cell
    coordinates = horzcat(X,Y); % save x and y as coordinates
    r_mat = squareform(pdist(coordinates)); % create pdist matrix for distances (r)
    
    gf = repmat(reshape(b,1,[]),size(b,1)^2,1); % assign gf values by recreating b matrix and repeating for each cell
    gf(gf > 0) = 1; % if value in gf is greater than 0 the value is assumed to be 1
    gf(gf <= 0) = 0; % if value in gf is equal to (or less than) 0 the value is assumed to be zero
    
    gc = repmat(w(:,t),[1,size(w(:,t))]);  % assign gc values from creating w for each column
    gc(gc > 0) = 1; % if there is water present (>0) the value is assuemd to be 1
    gc(gc <=0) = 0; % if there is no water present (<=0) the value is assumed to be 0
    
    g_max = ones(size(b,1)^2); % gmax matrix for both K and T representing maximum, corresponding to all 1s
    
    % Create f function to sum
    f_real_K = gf.*exp(-1*kf.^2.*(r_mat.^2)); % calculate f for every cell
    
    f_max_K = g_max.*exp(-1*kf.^2.*(r_mat.^2)); % calculate max f for every cell
    f_max_vec_K = sum(f_max_K,2); % sum max f across the rows
    f_max_mat_K = repmat(f_max_vec_K,1,size(b,1)^2); % replicate max f sum for every column
    
    f_K = f_real_K./f_max_mat_K; % normalize f real by f max
    sum_f_K = sum(f_K,2); % sum f across the rows
    
    % Calculate Ki from f and transform to original scene
    K = K0 + Kv.*(sum_f_K); 
    Kmat(:,t) = K;
    
    %%%%%% Ti 
    f_real_T = gc.*exp(-1*kc.^2.*(r_mat.^2)); % calculate f real for T
    
    f_max_T = g_max.*exp(-1*kc.^2.*(r_mat.^2)); % calculate f max for T
    f_max_vec_T = sum(f_max_T,2); % sum f max across rows
    f_max_mat_T = repmat(f_max_vec_T,1,size(b,1)^2); % replicate f max sum for every column
    
    f_T = f_real_T./f_max_mat_T; % normalize f real by f max
    sum_f_T = sum(f_T,1); % sum f across rows
    
    % Calculate Ti from f and transform to original scene
    T = Tmax.*(sum_f_T);
    Tmat(:,t) = T; % save for every time step
    %% Evaporation, Ei
    E = Ev.*bvec; % calculate E for each vegetated cell (bvec=1)
    indexE = find(bvec==0); % find bare cells (bvec = 0)
    E(indexE) = E(indexE) + Eb; % set E for each bare cell
    Emat(:,t) = E; % save E for each time step
    %% Biomass Growth
    indexkill = find(T<Tc); % find cells where T is less than threshold (Tc)
    b(indexkill) = b(indexkill)-growth; % decrease T<Tc cells by growth increment
    indexgrowth = find(T>=Tc); % find cells where T meets threshold (Tc)
    b(indexgrowth) = b(indexgrowth)+growth; % add T>=Tc cells by growth increment
    indexzero = find(b<0); % find cells where biomass is less than 0
    b(indexzero) = 0; % set negative biomass cells to 0
   
    wmat = reshape(w(:,t),size(b,1),size(b,2)); % reshape w into matrix in same dimensions as b
    
    indexw = find(wmat>wc & b==0); % find cells where water exceeds threshold AND there is no biomass
    b(indexw) = b(indexw) + growth; % add biomass growth increment to indexed cells
    wmat(indexw) = wmat(indexw) - wc; % subtract water threshold from water balance at indexed cells
    w(:,t) = reshape(wmat,1,numel(wmat)); % reshape water balance into vector 
    
    b(b > blimit) = blimit; % set biomass cells greater than limit to limit
    bmat(:,t) = reshape(b,1,[]); % save biomass as vector in matrix for each time step
    
    indexwater = find(w<0); % find negative values in the water balance for each step
    w(indexwater) = 0; % set negative water balance cells to zero
    
    % Sum water balance and recreate biomass scene
    w(:,t+1) = w(:,t) + (Imat(:,t) - Emat(:,t) - Tmat(:,t))*dt;
    
    % DO IT AGAIN FOR t+1????
    indexwater = find(w<0); % find negative values in the water balance for each step
    w(indexwater) = 0; % set negative water balance cells to zero
    
    % set probability values for next time step based on current time step
    p = min((K./P),1);
    
    
end

