% Function to call LCS Tool for calculating hyperbolic LCSs
% Can call this function to run FTLE and geodesics on the input data in
% either forward or backward time, with specified parameters

function [hFigure, ftle_, ftle_bkwd] = calculate_lcs_layered(data, domain, start_time, duration_T, duration_unit, resolution_x)

% calculate fwd or backward FTLE field and repelling or attracting hyperbolic LCSs using Haller's LCS tool
% data must follow Haller format conventions
% 
% INPUTS:
% domain must be given as [lon_min, lon_max ; lat_min, lat_max]
% start_time must be a datetime
% duration_T of time to integrate for
% duration_unit = 'hours' or 'days' --> specify what units the duration is in
% resolution_x: x-axis resolution for integration.
% dir: specify whether to calculate forward time, backward time, or both.
%       inputs allowed are "forward", "backward", and "both"
% 
% OUTPUTS:
% hFigure - figures of the FTLE field and the geodesic LCS
% optional additional outputs
% ftle_ - forward ftles calculated
% ftle_bkwd - backward ftles calculated


load(data);
date_dt = datetime(date(:, 2:7));
fignum = 1;
load('WEA_HFR_Data\wea_mboutline_lon_lat.mat')  % this line is specific to the WEA outline overlayed on the resulting graphs


% set up TIME parameters
start_ind = find(date_dt == start_time);
switch duration_unit
    case 'hours'
        end_time = start_time + hours(duration_T);
    case 'days'
        end_time = start_time + days(duration_T);
end
end_ind = find(date_dt == end_time);
assert(~isempty(end_ind));

start_id = time(start_ind); 
end_id = time(end_ind);
timespan = [start_id, end_id];
titledate = string(end_time, "yyyy-MM-dd, HH:00");

% set up DOMAIN parameters
resolutionX = resolution_x;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

% Parameter to skip Elliptic LCS calculations - this is always true for
% this version of the function, but needs to be there so it avoids sections
% of LCS Tool's functions that use this 
skip_elliptic = true;

% clear data not needed for the calculation to save memory space and remove issues with duplicate IDs for days
date = date(start_ind-1:end_ind+1,:);
time = time(start_ind-1:end_ind+1,:);
vLon = vLon(start_ind-1:end_ind+1, :, :);
vLat = vLat(start_ind-1:end_ind+1, :, :);

% % set up the gridded interpolant and specify the derivative function
% (derivative here is just calculating the gridded interpolant since we're using a real dataset and not a mathematically-generated flow)
interpMethod = 'spline';

vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);

lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
incompressible = true;  

%% LCS parameters

% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = .01; %default value - this should be between 0 and 0.5 and is recommended to be between 0.01 and 0.05 (i.e. 1-5% of grid spacing)

% Shrink lines (repelling LCSs)
shrinkLineMaxLength = 20;
shrinkLineLocalMaxDistance = 2*deltaX;
shrinkLineOdeSolverOptions = odeset('relTol',1e-6);

% Stretch lines (attracting LCSs)
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 4*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-6);

% Graphic properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];
initialPositionMarkerSize = 2;

% set up base plot (for Repelling and Attracting LCSs)
[hAxes, hFigure(fignum)] = setup_figure(domain, fignum);
title(hAxes, strcat('Repelling and Attracting LCSs:', {' '}, titledate))
subtitle(strcat(string(duration_T), {' '}, strrep(duration_unit, 's', ''), {' '}, 'integration', {', '}, 'Resolution:', {' '}, string(resolutionX), 'x', string(resolutionY)  ))
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')
fignum = fignum + 1;
   
%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

% Plot forward finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
%colormap(hAxes,flipud(gray))
colormap(hAxes,flipud(bone))
drawnow

%% Hyperbolic attracting LCSs
% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
[stretchLine,stretchLineInitialPosition] = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Plot hyperbolic attracting LCSs
hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
hAttractingLcs = [hAttractingLcs{:}];
set(hAttractingLcs,'color',attractingColor)
hStretchLineInitialPosition = arrayfun(@(idx)plot(hAxes,stretchLineInitialPosition(1,idx),stretchLineInitialPosition(2,idx)),1:size(stretchLineInitialPosition,2),'UniformOutput',false);
hStretchLineInitialPosition = [hStretchLineInitialPosition{:}];
set(hStretchLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
set(hStretchLineInitialPosition,'marker','o')
set(hStretchLineInitialPosition,'MarkerEdgeColor','w')
set(hStretchLineInitialPosition,'MarkerFaceColor',attractingColor)
drawnow
% hold on
% plot(wea_mboutline_lon_lat(:, 1), wea_mboutline_lon_lat(:, 2), 'k')

%% Hyperbolic repelling LCSs

% set up figure
% [hAxes, hFigure(fignum)] = setup_figure(domain, fignum);
% title(hAxes, strcat('Repelling LCSs:', {' '}, titledate))
% subtitle(strcat(string(duration_T), {' '}, strrep(duration_unit, 's', ''), {' '}, 'integration', {', '}, 'Resolution:', {' '}, string(resolutionX), 'x', string(resolutionY)  ))
% xlabel(hAxes,'Longitude (\circ)')
% ylabel(hAxes,'Latitude (\circ)')
% fignum = fignum + 1;

% Calculate repelling LCSs
[shrinkLine,shrinkLineInitialPosition] = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions);

% Plot hyperbolic repelling LCSs
hRepellingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),shrinkLine,'UniformOutput',false);
hRepellingLcs = [hRepellingLcs{:}];
set(hRepellingLcs,'color',repellingColor)
hShrinkLineInitialPosition = arrayfun(@(idx)plot(hAxes,shrinkLineInitialPosition(1,idx),shrinkLineInitialPosition(2,idx)),1:size(shrinkLineInitialPosition,2),'UniformOutput',false);
hShrinkLineInitialPosition = [hShrinkLineInitialPosition{:}];
set(hShrinkLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
set(hShrinkLineInitialPosition,'marker','o')
set(hShrinkLineInitialPosition,'MarkerEdgeColor','w')
set(hShrinkLineInitialPosition,'MarkerFaceColor',repellingColor)
drawnow
hold on
plot(wea_mboutline_lon_lat(:, 1), wea_mboutline_lon_lat(:, 2), 'k')
drawnow


% %% plot both
% % LH add - plot both together
% [hAxes, hFigure(fignum)] = setup_figure(domain, fignum);
% title(hAxes, strcat('Attracting and Repelling LCSs:', {' '}, titledate))
% subtitle(strcat(string(duration_T), {' '}, strrep(duration_unit, 's', ''), {' '}, 'integration', {', '}, 'Resolution:', {' '}, string(resolutionX), 'x', string(resolutionY)  ))
% xlabel(hAxes,'Longitude (\circ)')
% ylabel(hAxes,'Latitude (\circ)')
% fignum = fignum + 1;
% 
% % Plot forward finite-time Lyapunov exponent
% cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
% ftle_ = ftle(cgEigenvalue2,diff(timespan));
% plot_ftle(hAxes,domain,resolution,ftle_);
% %colormap(hAxes,flipud(gray))
% colormap(hAxes,flipud(bone))
% drawnow
% 
% % Copy objects from repelling LCS plot
% hRepellingLcs = copyobj(hRepellingLcs, hAxes);
% hAttractingLcs = copyobj(hAttractingLcs, hAxes);
% drawnow
% hold on
% line(wea_mboutline_lon_lat(:, 1), wea_mboutline_lon_lat(:, 2), 'color', 'k')
