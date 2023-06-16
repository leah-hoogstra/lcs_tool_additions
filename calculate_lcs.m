% Function to call LCS Tool for calculating hyperbolic LCSs
% Can call this function to run FTLE and geodesics on the input data in
% either forward or backward time, with specified parameters

function [hFigure, ftle_] = calculate_lcs(data, domain, start_time, duration, duration_unit, backward, resolution_x)

% calculate fwd or backward FTLE field and repelling or attracting hyperbolic LCSs using Haller's LCS tool
% data must follow Haller format conventions
% 
% INPUTS:
% domain must be given as [lon_min, lon_max ; lat_min, lat_max]
% start_time must be a datetime
% duration of time to integrate for
% duration_unit = 'hours' or 'days' --> specify what units the duration is in
% backward = true/false. default is false. use backward = true to integrate backward in time
% resolution: x-axis resolution for integration. 
% 
% OUTPUTS:
% hFigure - figures of the FTLE field and the geodesic LCS
% ftle_ - ftles calculated
% cg_Eigenvalue
% cg_Eigenvector

load(data);
date_dt = datetime(date(:, 2:7));
fignum = 1;

% set up TIME parameters
start_ind = find(date_dt == start_time);
switch duration_unit
    case 'hours'
        end_time = start_time + hours(duration);
    case 'days'
        end_time = start_time + days(duration);
end
end_ind = find(date_dt == end_time);
assert(~isempty(end_ind));

start_id = time(start_ind); 
end_id = time(end_ind);
timespan = [start_id, end_id];
if backward
   timespan = fliplr(timespan);  % flip times around if running in backward time
end

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

% set up base plot
if backward     % backward time --> attracting LCSs
    [hAxes, hFigure(fignum)] = setup_figure(domain, fignum);
    title(hAxes, 'Attracting LCSs')
    xlabel(hAxes,'Longitude (\circ)')
    ylabel(hAxes,'Latitude (\circ)')
    fignum = fignum + 1;
else            % forward time --> repelling LCSs
    [hAxes, hFigure(fignum)] = setup_figure(domain, fignum);
    title(hAxes, 'Repelling LCSs')
    xlabel(hAxes,'Longitude (\circ)')
    ylabel(hAxes,'Latitude (\circ)')
    fignum = fignum + 1;
end
   
%% Cauchy-Green strain eigenvalues and eigenvectors
if backward
    % Calculate CG Tensor eigenvals and eigenvecs for each point in domain
    % for backwards time, do this once with timespan backwards to get
    % backward in time FTLE field, then recalculate with time forwards
    % since Attracting LCS code is written to find them in forward time
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);
    
    % Plot finite-time Lyapunov exponent
    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    ftle_ = ftle(cgEigenvalue2,diff(timespan));
    plot_ftle(hAxes,domain,resolution,ftle_);
    %colormap(hAxes,flipud(gray))
    colormap(hAxes,flipud(bone))
    drawnow
    
    % for backwards time, recalculate using timespan flipped back to forward since the
    % geodesic Attracting LCS code is written to use forward time
    timespan = fliplr(timespan);
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);
else
    % Calculate CG Tensor eigenvals and eigenvecs for each point in domain
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);
    
    % Plot finite-time Lyapunov exponent
    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    ftle_ = ftle(cgEigenvalue2,diff(timespan));
    plot_ftle(hAxes,domain,resolution,ftle_);
    %colormap(hAxes,flipud(gray))
    colormap(hAxes,flipud(bone))
    drawnow
end



%% Hyperbolic repelling LCSs
if ~backward
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
end

%% Hyperbolic attracting LCSs
if backward
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
end