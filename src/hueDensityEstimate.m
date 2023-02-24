function pctHueTuned = hueDensityEstimate(rgcPerCone, sConeDensity)
% HUEDENSITYESTIMATE
%
% Description: 
%   Calculates the percentage of hue-tuned cells in the following way:
%   pctHueTuned = rgcPerCone * sConeDensity * numHueTunedTypes
%
% Syntax:
%   huePct = hueDensityEstimate(rgcPerCone, sConeDensity)
%
% Inputs:
%   rgcPerCone          double, scalar or vector
%       RGC:cone ratio
%   sConeDensity        double, scalar or vector (0-1)
%       Density of S-cones (S-cones/total cones)
% Output:
%   huePct              double, scalar, vector or matrix
%   The estimated percent of ganglion cells that are hue-tuned.
%   If a range of values are provided for one or both inputs, a range of
%   density estimates will be produced.
%
% Assumptions:
%   - Chromatic acuity is limited by S-cone density so each S-cone should
%     be represented by the full complement of hue-tuned cells. 
%   - We need 4 hue-tuned cells (M-LS, L-MS, LS-M and MS-L)

% History:
%   20Feb2023 - SSP
% -------------------------------------------------------------------------

    % If S-cone density wasn't provided as a decimal, make it one
    if sConeDensity > 1
        sConeDensity = sConeDensity/100;
    end

    % Let's say we have 100 cones
    numCones = 100;
    % And we want 4 hue-tuned cells per S-cone
    numHueTunedTypes = 4;


    % How many RGCs do we have?
    numRGCs = numCones * rgcPerCone;

    % How many S-cones do we have?
    numSCones = numCones .* sConeDensity;

    % How many hue-tuned cells do we have?
    numHueCells = numHueTunedTypes .* numSCones;

    % What percent of the RGCs are hue-tuned?
    pctHueTuned = 100 .* numHueCells' ./ numRGCs;

    % If a range of values were provided, print mean, SD, N and range
    if numel(pctHueTuned) > 3
        printStat(pctHueTuned(:), true);
        fprintf('Range = %.2f - %.2f\n', min(pctHueTuned(:)), max(pctHueTuned(:)));
    end