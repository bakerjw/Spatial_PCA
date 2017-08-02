function [sim_results] = simulateResiduals(x, y, T_sim, nsims, nPCs)

%% Simulation of within-event ground motion residuals at multiple periods using PCA
% Modified by Maryia Markhvida, 08/01/2017: added capability to simulate normalized 
%              within-event residuals using different number of principal components
% Created by Maryia Markhvida and Luis Ceferino, 09/13/2016
% Further methodology documentation: M. Markhvida, L. Ceferino, and J. Baker. “Modeling
%              spatially correlated spectral accelerations at multiple periods using 
%              principal component analysis and geostatistics”. Submitted to: Earthquake Engineering 
%               & Structural Dynamics In Review (2017).

% INPUT
% x, y     = x and y coordinates (in kilometers) for n locations of interest with size [n x 1]
% T_sim    = a vector specifying periods to be simulated of size [1 x t]
% nsims    = number of simulation to be performed
% nPCs     = number of principal components to be used in simulation
%           (recommended not less than 5)

% OUTPUT
% sim_results = [n x t] matrix of spatially cross-correlated within-event
%               normal residuals (mean 0 and variance 1) for n locations
%               and t periods

nLocs = length(x);

% Load PCA variogram model, PCA coefficients, and scaling factors
load variogramModel_19PC
load PCA_coefficients
load variance_scale_factor

% Create a matrix containing distances between locations
distanceMatrix = getDistanceMatrix(x,y);

% Scale variance if less than 19 principal components are used
if nPCs < 19
    for i = 1:length(modelVario)
        if strcmp(modelVario{i}.Type,'nug')
            modelVario{i}.Cn = modelVario{i}.Cn /variance_scale_factor(nPCs);
        else
            modelVario{i}.Cn = modelVario{i}.Cn /variance_scale_factor(nPCs);
            modelVario{i}.C1 = modelVario{i}.C1 /variance_scale_factor(nPCs);
            modelVario{i}.C2 = modelVario{i}.C2 /variance_scale_factor(nPCs);
        end
    end
end

% Create a covariance matricies for each of the principal components (PC's)
covMatrix = cell(1,nPCs);
for i = 1:nPCs
    covMatrix{i} =  getCovariance(distanceMatrix,modelVario{i});
end

% Simulate each of the PC's
sim_PCA = zeros(nLocs,nsims,nPCs);
mu = zeros(1,nLocs);
for i = 1:nPCs
    sim_PCA(:,:,i) = mvnrnd(mu,covMatrix{i},nsims)';
end

% Transform simulated PC's to spectral acceleration residuals for desired
% periods

[~,indexT] = ismember(T_sim,T);
sim_results = zeros(nLocs,length(indexT),nsims);


for i = 1:length(indexT)
    for j = 1:nsims
        temp_sim_PCA = reshape(sim_PCA(:,j,:),[nLocs,nPCs]);
        if indexT(i)
            sim_results(:,i,j) = PCA_T_X(temp_sim_PCA,PCAcoefs(indexT(i),1:nPCs),0);
        else
            for k = 1:nPCs
                extraPCAcoefs(1,k) = interp1(T',PCAcoefs(:,k),T_sim(i));
            end
            sim_results(:,i,j) = PCA_T_X(temp_sim_PCA,  extraPCAcoefs,0);
        end
    end
end

end

function [COV] = getCovariance(DIST,modelVario)

switch modelVario.Type
    case 'iso nest'
        [COV] = getIsoNestedCOV(modelVario,DIST);
    case 'nug'
        [COV] = getNugCOV(modelVario,DIST);
    otherwise
        disp('NOT RIGHT MODEL')
end
end