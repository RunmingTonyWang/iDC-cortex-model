%% Computational model of iDC‐induced neuronal polarization
% Author: Runming Wang, Gene Fridman
% 07/31/2025

%% Section 1: Plot Ve due to iDC and ΔVm along different neuron types (Fig.3AB)
clc; clear;

% 1) SETUP: GRID & ELECTRODE PARAMETERS
rho = 5e3;          % tissue resistivity (Ω·mm)
dx  = 0.02;        % horizontal resolution (mm)
dy  = 0.005;       % vertical resolution (mm)
CorticalDistance = -1.1:dx:1.1;  % horizontal range (mm)
Depth            =  0.0:dy:2.2;  % depth range (mm)
[X, Y] = meshgrid(CorticalDistance, Depth);

% iDC properties
Icenter_amp             = 20e-6;           % iDC current (A)
electrode_radius        = 0.125;           % iDC radius (mm)
electrode_center = [0, 0, 0];              % (x,y,z) center (mm)


% 2) EXTRACELLULAR POTENTIAL (Ve)
Ve  = zeros(numel(Depth), numel(CorticalDistance));
els = 0;
for theta = 1:360
    for dist = 0:0.01:electrode_radius
        x0 = cosd(theta)*dist + electrode_center(1);
        y0 =                  electrode_center(2);
        z0 = sind(theta)*dist + electrode_center(3);
        r  = sqrt((X - x0).^2 + (Y - y0).^2 + z0.^2) + 1e-5;
        Va = rho ./ (2*pi*r);
        Ve = Ve + Va;
        els = els + 1;
    end
end
Ve = Ve * Icenter_amp / els; % Scale by total current


% 3) DEFINE NEURON TYPES
% format: name, topDepth (mm), bottomDepth (mm), somaDepth (mm)
neuronTypes = [
    % Excitatory neurons
    struct('name','L2py',       'top',0.03, 'bottom',0.60, 'soma',0.27)
    struct('name','L3py',       'top',0.03, 'bottom',1.05, 'soma',0.44)
    struct('name','L4sp',       'top',0.20, 'bottom',1.10, 'soma',0.68)
    struct('name','L4ss',       'top',0.40, 'bottom',1.00, 'soma',0.73)
    struct('name','L5st',       'top',0.03, 'bottom',1.45, 'soma',1.14)
    struct('name','L5tt',       'top',0.03, 'bottom',1.45, 'soma',1.21)
    struct('name','L6cc',       'top',0.75, 'bottom',1.75, 'soma',1.43)
    struct('name','L6ct',       'top',0.70, 'bottom',2.10, 'soma',1.53)
    % Inhibitory neurons
    struct('name','BasketCell1','top',0.55, 'bottom',0.75, 'soma',0.65)
    struct('name','BasketCell2','top',1.40, 'bottom',1.90, 'soma',1.65)
    struct('name','Martinotti1','top',0.25, 'bottom',0.55, 'soma',0.40)
    struct('name','Martinotti2','top',1.00, 'bottom',1.75, 'soma',1.20)
];


% 4) MIRROR ESTIMATE FOR EACH NEURON TYPE
PNm_all = cell(numel(neuronTypes),1);
for i = 1:numel(neuronTypes)
    PNm_this = zeros(size(Ve));
    topIdx    = find(Depth >= neuronTypes(i).top, 1);
    bottomIdx = find(Depth <= neuronTypes(i).bottom, 1, 'last');
    for col = 1:3:numel(CorticalDistance) % using a step of 3 for visualization
        segment    = Ve(topIdx:bottomIdx, col);
        meanSegment = mean(segment);
        PNm_this(topIdx:bottomIdx, col) = meanSegment - segment;
    end
    PNm_all{i} = PNm_this;
end


% 5) VISUALIZATION: VE AND PNM
Ve_lims  = [-0.03, 0.03];           % For Ve (V)
PNm_lims = [-0.015, 0.015];         % For Vm changes (V)
layers      = [0.15, 0.3, 0.55, 0.9, 1.4, 1.95];    %layer boundaries (mm)
%layerLabels = {'II','III','IV','V','VI','WM'};      %layer labels

cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
        ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];

figure('Color','w');
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');

% Ve subplot (Fig.3A)
ax=nexttile;
imagesc(CorticalDistance,Depth,Ve,Ve_lims); colormap(ax,'parula'); caxis(Ve_lims);
colorbar; title('V_e (V)'); xlabel('Lateral (mm)'); ylabel('Depth (mm)'); axis equal tight; 
set(ax,'TickDir','out'); hold on;
for L=layers; plot([min(CorticalDistance) max(CorticalDistance)], [L L],'--k'); end; hold off;

% PNm subplots (Fig.3B)
for i=1:numel(neuronTypes)
    ax=nexttile;
    imagesc(CorticalDistance,Depth,PNm_all{i},PNm_lims); colormap(ax,cmap); caxis(PNm_lims);
    colorbar; title(neuronTypes(i).name,'Interpreter','none');
    xlabel('Lateral (mm)'); ylabel('Depth (mm)'); axis equal tight; set(ax,'TickDir','out'); hold on;
    for L=layers; plot([min(CorticalDistance) max(CorticalDistance)], [L L],'--k'); end;
    % Soma markers
    somaCols = 1:3:numel(CorticalDistance);
    plot(CorticalDistance(somaCols), neuronTypes(i).soma*ones(size(somaCols)),'ko','LineWidth',1.5);
    hold off;
end




%% Section 2: Weighted 2D Interpolation and Layer-wise Average Versus Distance (Fig.3CDE)
clc; clear;

% 1) SETUP: GRID & ELECTRODE PARAMETERS
rho = 5e3;        % tissue resistivity (Ω·mm)
dx = 0.001;       % horizontal resolution (mm)
dy = 0.005;       % vertical resolution (mm)
CorticalDistance = -2:dx:2;  % horizontal range (mm)
Depth = 0.0:dy:2.2;          % depth range (mm)
[X, Y] = meshgrid(CorticalDistance, Depth);

% iDC properties
Icenter_amp = +20e-6;                    % iDC current (A)
electrode_radius = 0.125;         % iDC radius (mm)
electrode_center = [0, 0, 0];     % (x,y,z) center (mm)


% 2) CALCULATE EXTRACELLULAR POTENTIAL Ve
Ve = zeros(length(Depth), length(CorticalDistance));
els = 0;
for theta = 1:360
    for dist = 0:0.01:electrode_radius
        x0 = cosd(theta)*dist + electrode_center(1);
        y0 = electrode_center(2);
        z0 = sind(theta)*dist + electrode_center(3);
        distances = sqrt((X - x0).^2 + (Y - y0).^2 + z0.^2) + 1e-4;  % avoid division by zero
        Va = rho ./ (2*pi*distances);
        Ve = Ve + Va;
        els = els + 1;
    end
end
Ve = Ve * Icenter_amp / els; % Scale by total current


% 3) DEFINE NEURON TYPES & RELATIVE DENSITIES
% format: name, topDepth (mm), bottomDepth (mm), somaDepth (mm)
neuronTypes = [
    struct('name','L2py', 'topDepth',0.03,  'bottomDepth',0.60, 'somaDepth',0.27)
    struct('name','L3py', 'topDepth',0.03,  'bottomDepth',1.05, 'somaDepth',0.44)
    struct('name','L4sp', 'topDepth',0.20, 'bottomDepth',1.10, 'somaDepth',0.68)
    struct('name','L4ss', 'topDepth',0.40, 'bottomDepth',1.00, 'somaDepth',0.73)
    struct('name','L5st', 'topDepth',0.03,  'bottomDepth',1.45, 'somaDepth',1.14)
    struct('name','L5tt', 'topDepth',0.03,  'bottomDepth',1.45, 'somaDepth',1.21)
    struct('name','L6cc', 'topDepth',0.75, 'bottomDepth',1.75, 'somaDepth',1.43)
    struct('name','L6ct', 'topDepth',0.70, 'bottomDepth',2.10, 'somaDepth',1.53)
];
% Relative density percentages for each neuron type
densities_percent = [7.51, 11.01, 17.70, 18.93, 6.89, 11.11, 11.52, 15.33]; 
densities_fraction = densities_percent / sum(densities_percent);
% Total number of rods to place (across all neuron types)
N_rods_total = 4000;
rods_per_type = round(densities_fraction * N_rods_total);


% 4) RANDOM SCATTERING NEURON RODS WITHOUT OVERLAP 
PNm_random = zeros(size(Ve));
possibleCols = 1:length(CorticalDistance);
usedCols = false(size(possibleCols));  % track used columns
somaX = [];
somaY = [];
for iType = 1:length(neuronTypes)
    topIdx    = find(Depth >= neuronTypes(iType).topDepth, 1, 'first');
    bottomIdx = find(Depth <= neuronTypes(iType).bottomDepth, 1, 'last');
    if neuronTypes(iType).somaDepth < neuronTypes(iType).topDepth || ...
       neuronTypes(iType).somaDepth > neuronTypes(iType).bottomDepth
       warning('Soma depth for %s is outside the rod range!', neuronTypes(iType).name);
    end
    for r = 1:rods_per_type(iType)
        picked = false;
        while ~picked
            randIdx = randi([1, length(possibleCols)]);
            if ~usedCols(randIdx)
                usedCols(randIdx) = true;
                picked = true;
            end
        end
        xIdx = possibleCols(randIdx);
        % Compute mirror estimate for this rod (vertical segment)
        Ve_segment = Ve(topIdx:bottomIdx, xIdx);
        Ve_mean = mean(Ve_segment);
        PNm_segment = Ve_mean - Ve_segment;
        PNm_random(topIdx:bottomIdx, xIdx) = PNm_segment;
        % Record soma location
        somaX(end+1) = CorticalDistance(xIdx);
        somaY(end+1) = neuronTypes(iType).somaDepth;
    end
end
PNm_lims = [-0.015, 0.015];
figure('Name','Randomly Scattered Neurons Combined','Color','w');
imagesc(CorticalDistance, Depth, PNm_random, PNm_lims);
axis equal tight;
set(gca, 'TickDir','out');
colorbar;
title('Combined Randomly Scattered Neurons');
xlabel('Horizontal Distance (mm)');
ylabel('Depth (mm)');
cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
        ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
colormap(gca, cmap);
caxis(PNm_lims);
layers = [0.15, 0.3, 0.55, 0.9, 1.4, 1.95];  %Overlay cortical layer boundaries
% layer_labels = {'II','III','IV','V','VI', 'WM'};
hold on;
for L = 1:length(layers)
    plot([min(CorticalDistance), max(CorticalDistance)], [layers(L) layers(L)], '--k');
end
plot(somaX, somaY, 'ko', 'MarkerSize',6, 'LineWidth',1.5); % Plot soma markers
hold off;


% 5) WEIGHTED-SUM INTERPOLATION OF SOMA POTENTIALS (Fig.3CD)
nNeurons = length(somaX);
somaPotential = zeros(1, nNeurons);
for i = 1:nNeurons
    [~, colIdx] = min(abs(CorticalDistance - somaX(i)));
    [~, rowIdx] = min(abs(Depth - somaY(i)));
    somaPotential(i) = PNm_random(rowIdx, colIdx);
end
fineX = linspace(min(CorticalDistance), max(CorticalDistance), 200);
fineY = linspace(min(Depth), max(Depth), 200);
[XX, YY] = meshgrid(fineX, fineY);
sigma = 0.1;  % Kernel bandwidth in mm
weightedSum = zeros(size(XX));
for i = 1:length(somaX)
    dist2 = (XX - somaX(i)).^2 + (YY - somaY(i)).^2;
    w = exp(-dist2/(2*sigma^2));
    weightedSum = weightedSum + w * somaPotential(i);
end
interpSum = weightedSum; 

figure('Name','2D Weighted Sum of Soma Potentials','Color','w');
imagesc(fineX, fineY, interpSum);
set(gca, 'YDir', 'reverse');
set(gca, 'TickDir','out');
axis equal tight;
caxis(PNm_lims * 50); %For visualization
colormap(cmap);
colorbar;
xlabel('Horizontal Distance (mm)');
ylabel('Cortical Depth (mm)');
yticks(0:0.1:2.2);
ylim([0 2.2]);
title('2D Weighted Sum of Soma Potentials');
hold on;
for L = 1:length(layers)
    plot([min(fineX) max(fineX)], [layers(L) layers(L)], 'k--', 'LineWidth', 1);
end
hold off;


% 6) LAYER-WISE AVERAGE VS DISTANCE (Fig.3E)
% Distances of interest (mm)
hDists = [0.2, 0.55, 0.9, 1.25];
nDist  = numel(hDists);
% Layer depth ranges (mm)
layerRanges = [
    0.15, 0.55;   % Layer 2/3
    0.55, 0.90;   % Layer 4
    0.90, 1.40;   % Layer 5
    1.40, 1.65    % Layer 6
];
layerNames = {'L2/3','L4','L5','L6'};
nLayers    = size(layerRanges,1);

avgPot = zeros(nLayers, nDist);
xIdx = arrayfun(@(d) find(abs(fineX - d)==min(abs(fineX-d)),1), hDists);
for L = 1:nLayers
    yIdx = find(fineY >= layerRanges(L,1) & fineY <= layerRanges(L,2));
    for k = 1:nDist
        avgPot(L,k) = mean( interpSum(yIdx, xIdx(k)) );
    end
end

figure('Name','Layer-wise Soma Potential vs Distance','Color','w');
hold on;
colors = lines(nLayers);
for L = 1:nLayers
    plot(hDists, avgPot(L,:), '-o', ...
         'Color', colors(L,:), ...
         'LineWidth', 2, ...
         'DisplayName', layerNames{L});
end
hold off;
xlabel('Horizontal Distance from Catheter (mm)');
ylabel('Average Weighted Soma Potential (V)');
legend('Location','best');
title('Depth-Layer Averages of Weighted Soma Potentials');
grid on;
ylim([-0.5 0.5]);
xticks(hDists);
