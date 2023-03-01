function rgcPerCone = rgcPerConeEstimate(sConeDensity, pctMidgetRGC)
% RGCPERCONEESTIMATE
%
% Syntax:
%   rgcPerCone = rgcPerConeEstimate(sConeDensity, pctMidgetRGC)
%
% Assumptions: 
% - The percentage of foveal RGCs that are midget RGCs is ~86%. This comes
%   from transcriptomics work by Peng et al (2019). Given the size of their
%   dataset and the minimal possibility of sampling bias in their methods
%   (see paper for discussion on this), their estimate seems reliable. 
%   Additional confidence comes from our own unpublished work where we find
%   that at least 16% of foveal RGCs had response properties inconsistent 
%   with midget RGCs. 
% -------------------------------------------------------------------------

    if nargin < 2
        pctMidgetRGC = 0.86;    % Peng et al., 2019 (see above)
    end

    % Let's say we have 1000 cones, starting with a # makes the logic clear
    numCones = 1000;

    % This population can be divided into L/M-cones and S-cones
    numConesLM = (1-sConeDensity) * numCones;
    numConesS = sConeDensity * numCones;

    % We know there are at least 2 midget RGCs per L/M-cone (OFF and ON)
    numMidgetRGCs = 2 * numConesLM;
    % And 1 midget RGC per S-cone (OFF)
    numMidgetRGCs = numMidgetRGCs + numConesS;

    % Then we should now be able to obtain the total number of RGCs. 
    numRGCs = numMidgetRGCs' ./ pctMidgetRGC;

    % Divide by the number of cones and we have our ratio
    rgcPerCone = numRGCs / numCones;
    



    