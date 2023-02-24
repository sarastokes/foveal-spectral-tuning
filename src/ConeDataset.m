classdef (Abstract) ConeDataset < handle 
% CONEDATASET (Abstract)
%
% Description:
%   Abstract class encapsulating classification of ROIs based on their 
%   responses to sinusoidally-modulated chromatic stimuli. 
% 
% ROIs can be omitted for three reasons:
% - Bad segmentation (badSegIDs): these ROIs are not included in the total
%   number of cells analyzed as they are not valid cells
% - Unresponsive ('NR' in responseMap): these cells did not have an average 
%   d-prime over 2 SDs are are considered "unresponsive"
% - Unreliable (unreliableIDs): these ROIs were responsive but did not 
%   exceed 2 SDs in at least 2 experiments. There are many reasons why this 
%   could occur, but because one is crosstalk, they are not classified
% - isBad indexes ROIs omitted for any reason
%
% Ultimately ROIs need to be significantly responsive and reliable
% - passTable declares which ROIs meet these standards
% - stimResponsivity looks only at responsivity on a per-stimulus basis
% - stimReliability looks only at reliability on a per-stimulus basis
% - isUnresponsive indexes cells unresponsive to all stimuli
%
% ROIs are passed through several levels of classification:
% - responseMap identifies ROIs by stimulus responsivity (not exclusive)
% - comparisonMap identifies ROIs by 2 cone comparisons (not exclusive)
% - coneMap identifies ROIs by full cone combinations (exclusive)
%
% The final classification looks at cone opponency and polarity
% - roiClasses is a string array with final classes
% - achromClasses contains only the achromatic responses which are 
%   overwritten by cone opponency in roiClasses
%
% References:
%   Based on Tyler Godat's analysis scripts. Input data comes from Tyler's 
%   datasets and markup (e.g., unreliableIDs, badSegIDs)
%
% History:
%   9Feb2023 - SSP
% -------------------------------------------------------------------------

    properties (SetAccess = private)
        % Animal number 
        animalID            {mustBeInteger}
        % Tyler's proccessed dataset (details in constant props)
        dataset
        % These were unreliable (not significant in at least 2 trials)
        unreliableIDs 
        % These were improperly segmented and shouldn't be in total count
        badSegIDs
    end

    properties (SetAccess = private)
        % Unique identifier for each ROI (animalID.roiID)
        roiIDs
        % The total number of ROIs in the dataset
        numRois
        % The number of analyzed ROIs (all - bad)
        numAnalyzed 
        % The number of significantly responsive ROIs (all - bad - NR)
        numResponsive
        % The number of reliably responsive ROIs (all - bad - unreliable) 
        numReliable 
    end

    properties (SetAccess = private)
        % Average SNR for each ROI (numRois x numStim)
        avgSNR
        % Average phase for each ROI (numRois x numStim)
        avgPhase
        % Polarity for each ROI (numRois x numStim)
        avgSign
        % Average SNR with polarity sign (numRois x numStim)
        signedSNR

        stimResponsivity
        % Matrix designating reliable responsivity (numRois x numStim)
        stimReliability

    end

    properties (SetAccess = protected)
        % Table of response validity (numROIs x numStim)
        passTable
        % Cells that shouldn't be classified for any reason
        isBad
        % Cells that were not responsive (excluding bad cells)
        isUnresponsive

        % Responsivity combinations, no phase considered, exclusive
        responseMap
        % Cone comparisons, not exclusive
        comparisonMap
        % Classification by cone combinations, exclusive
        coneMap
        
        % Final classification
        roiClasses
        % Classification on achromatic stimulus alone
        achromClasses
    end

    properties (Hidden, Dependent)
        % IDs of omitted cells (derived from isBad)
        omittedIDs
        % Number of stimuli presented
        numStim
        % Number of experiments
        numReps 
    end

    properties (Hidden, Constant)
        DIMS = ["RoiIDs", "Stats", "Stimuli", "Repeats"];
        STATS = ["Phase", "Amplitude", "SNR", "Pass"];
    end

    properties (Abstract)
        STIMULI
    end

    methods 
        function obj = ConeDataset(animalID, dataset, badSegIDs, unreliableIDs)

            % Input assignment
            obj.animalID = animalID;
            obj.dataset = dataset;
            if nargin > 2 && ~isempty(badSegIDs)
                obj.badSegIDs = sort(badSegIDs);
            end
            if nargin > 3 && ~isempty(unreliableIDs)
                obj.unreliableIDs = sort(unreliableIDs);
            end

            % Determine sample characteristics
            obj.numRois = size(obj.dataset, 1);
            obj.numAnalyzed = obj.numRois - numel(obj.badSegIDs);
            fprintf('Analyzing %u of %u\n', obj.numAnalyzed, obj.numRois);

            % Assign unique IDs, integer=animalID, decimals=roiID (3 sig)
            obj.roiIDs = 1:obj.numRois;
            obj.roiIDs = (0.001 * obj.roiIDs) + obj.animalID;

            obj.analyzeMetrics();
            obj.classify();
        end
    end

    % Dependent set/get methods
    methods 
        function out = get.numReps(obj)
            if isempty(obj.dataset)
                out = 0;
            else
                out = size(obj.dataset, 4);
            end
        end

        function out = get.numStim(obj)
            if isempty(obj.dataset)
                out = 0;
            else
                out = size(obj.dataset, 3);
            end
        end

        function out = get.omittedIDs(obj)
            if isempty(obj.unreliableIDs) && isempty(obj.badSegIDs)
                out = [];
            else
                out = sort(unique([obj.unreliableIDs, obj.badSegIDs]));
            end
        end
    end

    % Base methods
    methods (Sealed)      
        function analyzeMetrics(obj)
            % Extract average SNR and phase
            obj.avgSNR = squeeze(mean(obj.dataset(:, 3, :, :), 4));
            obj.avgPhase = squeeze(rad2deg(circ_mean(...
                deg2rad(obj.dataset(:,1,:,:)-180), [], 4))+180);

            % Determine classify responses as ON- or OFF-dominant
            obj.avgSign = obj.avgPhase;
            obj.avgSign(obj.avgPhase < 180) = -1;
            obj.avgSign(obj.avgPhase >= 180) = 1;

            % SNR with + or - for ON/OFF
            obj.signedSNR = obj.avgSign .* obj.avgSNR;
        end

        function classify(obj)
            obj.reset();

            % Looks at which ROIs responded to each stimulus
            obj.calcStimResponsivity();

            % Groups ROIs by which cones had significant responses
            obj.calcConeCombinations();

            % Count up cells by attribute (responsive, reliable, etc)
            obj.calcCounts();

            % Non-exclusive classification of cone comparisons
            obj.calcComparisons();

            % Classifies different forms of opponency/non-opponency
            obj.calcOpponency();

            % Mark cells that weren't classified with reason for omission
            obj.classifyBad();

            % Further subdivides cone groupings by cone polarity
            obj.classifyPolarity();


            % Report results to the command line
            obj.reportResponsivity();
            obj.reportReliability();
            obj.reportResponseClasses();
            obj.reportOpponency();
            obj.reportPolarity();
        end
         
        function reset(obj)
            % Preallocation of classification arrays
            obj.roiClasses = repmat("", [obj.numRois, 1]);
            obj.achromClasses = repmat("", [obj.numRois, 1]);
            
            % Containers
            obj.responseMap = containers.Map();
            obj.comparisonMap = containers.Map();
            obj.coneMap = containers.Map();

            % Unset
            obj.passTable = [];
        end
    end

    methods
        function xy = getSACoords(obj, coneType)
            if nargin < 2
                roiID = 1:obj.numRois;
                roiID(obj.omittedIDs) = [];
            else              
                if isnumeric(coneType)
                    roiID = coneType;
                else
                    roiID = obj.coneMap(coneType);
                end
            end

            if isempty(roiID)
                xy = [];
                return
            end

            x = obj.signedSNR(roiID, obj.STIMULI == "A");
            y = obj.signedSNR(roiID, obj.STIMULI == "S");
            xy = [x, y];
        end
    end

    % Should be extended by subclasses
    methods (Access = protected)
        function calcConeCombinations(obj)
            obj.responseMap('NoResp')  = find(obj.isUnresponsive);
            obj.responseMap('NoResp') = setxor(obj.responseMap('NoResp'), obj.badSegIDs);
            obj.responseMap('Unrely') = obj.unreliableIDs;
            obj.responseMap('BadSeg') = obj.badSegIDs;
        end

        function calcOpponency(obj)
            error('Not yet implemented! Must be defined by subclasses');
        end

        function calcComparisons(obj)
            error('Not yet implemented! Must be defined by subclasses');
        end
    end

    methods
        function calcBadCells(obj)
            obj.isUnresponsive = sum(obj.passTable{:, obj.STIMULI ~= "C"}) == 0;
            obj.isUnresponsive(obj.badSegIDs) = 0;

            obj.isBad = obj.isUnresponsive;
            obj.isBad(obj.omittedIDs) = 1;
        end

        function calcStimResponsivity(obj)
            
            obj.passTable = array2table(obj.avgSNR >= 2,... 
                'VariableNames', obj.STIMULI);
            
            % Identify unresponsive cells
            obj.isUnresponsive = sum(obj.passTable{:, obj.STIMULI ~= "C"}, 2) == 0;
            % Create an index specifying all bad cells
            obj.isBad = obj.isUnresponsive;
            obj.isBad(obj.omittedIDs) = 1;
            
            % Zero unresponsive cells
            obj.passTable{obj.badSegIDs, :} = 0;
            obj.stimResponsivity = sum(obj.passTable{:,:}, 1);

            % Zero unreliable cells
            obj.passTable{obj.unreliableIDs, :} = 0;
            obj.stimReliability = sum(obj.passTable{:,:}, 1);
        end

        function calcCounts(obj)
            % Analyzed ROIs that responded significantly 
            obj.numResponsive = obj.numAnalyzed - numel(obj.responseMap('NoResp'));
            % Analyzed ROIs that responded significantly AND reliably 
            obj.numReliable = obj.numAnalyzed - numel(unique([obj.responseMap('NoResp'); obj.unreliableIDs']));
        end
    end

    % Classification methods
    methods 
        function classifyBad(obj)
            % Cells that register as significant but only because one large
            % response is skewing the average SNR for a stimulus
            obj.roiClasses(obj.unreliableIDs) = "Unrely";
            % These cells did not respond significantly to any stimulus
            obj.roiClasses(obj.responseMap('NoResp')) = "NR";
            % These cells were badly segmented and should be completely 
            % omitted from the analysis
            obj.roiClasses(obj.badSegIDs) = "OMIT";
        end

        function classifyPolarity(obj)
            % Start with achromatic, if present
            iStim = find(obj.STIMULI == "A");
            if ~isempty(iStim)
                iRois = obj.responseMap('A');
                for i = 1:numel(iRois)
                    % Check to see if ROI is omitted
                    if obj.roiClasses(iRois(i)) ~= ""
                        continue
                    end
                    if obj.avgSign(iRois(i), iStim) > 0
                        obj.roiClasses(iRois(i)) = "A-ON";
                    else
                        obj.roiClasses(iRois(i)) = "A-OFF";
                    end
                end
                % Save out as cone opponency below takes precedence
                obj.achromClasses = obj.roiClasses;
                % Redefine some classes for easy interpretation
                obj.achromClasses(ismember(obj.achromClasses, ["NR", "OMIT", "Unrely"])) = "Omit";
                obj.achromClasses(obj.achromClasses == "") = "NR";
            end

            % Continue with cone-opponent, overwriting achrom if needed
            k = obj.coneMap.keys;

            % Order by 1. number of cones, 2. opponency
            isOpponent = cellfun(@(x) contains(x, 'v'), k);
            [isOpponent, idx] = sort(isOpponent);
            k = k(idx);
            
            nCones = cellfun(@numel, k);
            nCones(isOpponent) = nCones(isOpponent) - 1;
            [~, idx] = sort(nCones);
            k = k(idx);

            % Assign polarity
            for i = 1:numel(k)
                iKey = k{i};

                % Find the stimulus index
                iStim = find(obj.STIMULI == string(iKey(1)));
                if isempty(iStim) && startsWith(iKey, ["LvM", "MvL"])
                    iStim = find(obj.STIMULI == "I");
                end 

                % Find the ROIs matching the key, if they exist
                iRois = obj.coneMap(iKey);
                if isempty(iRois)
                    continue
                end

                % Split into ON and OFF based on 1st cone response 
                % e.g. for LvMS, polarity is based on L-cone response
                for j = 1:numel(iRois)
                    
                    % Check to see if ROI is omitted
                    if obj.roiClasses(iRois(j)) ~= ""
                        continue
                    end

                    if obj.avgSign(iRois(j), iStim) > 0
                        obj.roiClasses(iRois(j)) = sprintf("%s-ON", iKey);
                        %fprintf('Set %u as %s\n', iRois(j), sprintf("%s-ON", iKey));
                    else
                        obj.roiClasses(iRois(j)) = sprintf("%s-OFF", iKey);
                        %fprintf('Set %u as %s\n', iRois(j), sprintf("%s-OFF", iKey));
                    end
                end
            end
        end
    end

    % Reporting methods
    methods 
        function reportResponsivity(obj)
            fprintf('\nTOTAL ANALYZED: %u of %u\n', obj.numAnalyzed, obj.numRois);
            fprintf('\nRESPONSIVITY\n');
            fprintf('Responsive cells %u of %u (%.2f%%).',...,...
                obj.numResponsive, obj.numAnalyzed,...
                100 * obj.numResponsive/obj.numAnalyzed);
            fprintf(' By stimulus:\n');
            for i = 1:numel(obj.STIMULI)
                fprintf('\t%s = %u (%.2f%%)\n', ...
                    obj.STIMULI(i), obj.stimResponsivity(i),...
                    obj.stimResponsivity(i)/obj.numAnalyzed * 100);
            end
            fprintf('\n');
        end

        function reportReliability(obj)
            fprintf('\nRELIABILITY\n');
            fprintf('Reliable cells %u of %u (%.2f%%).',...
                obj.numReliable, obj.numAnalyzed, ...
                100 * obj.numReliable/obj.numAnalyzed);
            fprintf(' By stimulus:\n');
            for i = 1:numel(obj.STIMULI)
                fprintf('\t%s = %u (%.2f%%)\n',... 
                    obj.STIMULI(i), obj.stimReliability(i),...
                    100 * obj.stimReliability(i)/obj.numAnalyzed);
            end
            fprintf('\n');
        end

        function reportOpponency(obj)
            
            % Extract opponency type from keys
            k = obj.coneMap.keys;

            % Find the opponent keys and sort by number of cones
            fprintf('\nOPPONENT\n');
            isOpponent = cellfun(@(x) contains(x, 'v'), k);
            isDichromatic = cellfun(@(x) numel(x) == 3, k(isOpponent));
            [~, kOrder] = sort(isDichromatic, 'descend');
            kIdx = find(isOpponent); kIdx = kIdx(kOrder);

            for i = 1:numel(kIdx)
                iName = k{kIdx(i)};
                fprintf('\t%s = %u\n', iName, numel(obj.coneMap(iName)));
            end

            % Find the non-opponent keys and sort by number of cones
            fprintf('\nNONOPPONENT\n');
            numCones = cellfun(@(x) numel(x), k(~isOpponent));
            [~, kOrder] = sort(numCones, 'descend');
            kIdx = find(~isOpponent); kIdx = kIdx(kOrder);

            for i = 1:numel(kIdx)
                iName = k{kIdx(i)};
                fprintf('\t%s = %u\n', iName, numel(obj.coneMap(iName)));
            end
        end

        function reportResponseClasses(obj)
            fprintf('\nRESPONSE CLASSIFICATION\n');

            k = obj.responseMap.keys;
            idx = cellfun(@(x) ismember(string(x), ["NoResp", "Unrely", "BadSeg"]), k);

            for i = 1:numel(k)
                if ~idx(i)
                    fprintf('\t%s = %u\n', k{i}, numel(obj.responseMap(k{i})));
                end
            end
            fprintf('\tNonresponsive = %u\n', numel(obj.responseMap('NoResp')));
            fprintf('\tUnreliable = %u\n', numel(obj.responseMap('Unrely')));
            fprintf('\tBad ROI = %u\n', nnz(obj.responseMap('BadSeg')));
            fprintf('\n');
        end

        function reportPolarity(obj)
            fprintf('\nPOLARITY\n');

            % All stimuli, except for the control stimulus
            stimLetters = obj.STIMULI(obj.STIMULI ~= "C");
            for i = 1:numel(stimLetters)
                fprintf('\t%s-ON = %u\n', stimLetters(i),...
                    nnz(obj.roiClasses == sprintf("%s-ON", stimLetters(i))));
                fprintf('\t%s-OFF = %u\n', stimLetters(i),...
                    nnz(obj.roiClasses == sprintf("%s-OFF", stimLetters(i))));
            end
            fprintf('\n');
        end
    end
end 