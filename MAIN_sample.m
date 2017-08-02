% sample main file for simulating within-event normal residuals
clear all; close all; clc;
addpath('./Functions')
addpath('./Input data')

T_sim = [0.01,0.1,1]; % vector of periods of interest in range 0.01-5s
nsims = 1000;% number of simulations
nPCs = 18;   % number of principal components to be used in simulation
             % Use nPCs = 18 (or 19) if comupational efficiency is sufficient 
             % Do not use nPC < 5;

x_coordinates = randi(200,100,1); % Should be [n x 1] in kilometers
y_coordinates = randi(200,100,1); % Should be [n x 1] in kilometers

[sim_results] = simulateResiduals(x_coordinates, y_coordinates, T_sim, nsims, nPCs);
