classdef TylerConeDataset < handle 

    properties (SetAccess = private)
        animalID 
        % Animal ID
        source
        % Tyler's proccessed dataset (details in constant props)
        dataset
        % These were unreliable (not significant in at least 2 trials)
        unreliableIDs 
        % These were improperly segmented and shouldn't be in total count
        badSegIDs
    end

    properties 
        numRois
        numAnalyzed 
        numReps
        roiIDs
    end

    properties (SetAccess = private)
        avgSNR
        avgPhase
        avgSign
        signedSNR

        stimResponsivity
        stimReliability

        roiClasses
        roiCones
        roiAchrom
    end

    properties (SetAccess = private)
        passTable
        isBad
        isUnresponsive
        numResponsive
        numReliable

        % Responsivity combinations, no phase considered, exclusive
        responseMap
        % Cone comparisons, groups are not exclusive
        comparisonMap
        % Classification by cone combinations, exclusive
        coneMap
    end

    properties (Hidden, Dependent)
        omittedIDs
    end

    properties (Hidden, Constant)
        DIMS = ["RoiIDs", "Stats", "Stimuli", "Repeats"];
        % Control, L, M, S, Achromatic
        STIMULI = ["C", "L", "M", "S", "A"];
        STATS = ["Phase", "Amplitude", "SNR", "Pass"];
    end

    methods 
        function obj = TylerConeDataset(source, animalID, dataset, badSegIDs, unreliableIDs)
            obj.source = source;
            obj.animalID = animalID;

            obj.dataset = dataset;
            obj.badSegIDs = sort(badSegIDs);
            obj.unreliableIDs = sort(unreliableIDs);

            % Determine sample characteristics
            obj.numReps = size(obj.dataset, 4);
            obj.numRois = size(obj.dataset, 1);
            obj.numAnalyzed = obj.numRois - numel(obj.badSegIDs);
            fprintf('Analyzing %u of %u\n', obj.numAnalyzed, obj.numRois);

            % Integer = animalID, decimals = cell ID (3 sig digits)
            obj.roiIDs = 1:obj.numRois;
            obj.roiIDs = (0.001 * obj.roiIDs) + obj.animalID;

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

            % Preallocation of classification arrays
            obj.roiClasses = repmat("", [obj.numRois, 1]);
            obj.roiCones = repmat("", [obj.numRois, 1]);
            obj.roiAchrom = repmat("", [obj.numRois, 1]);
        end

        function out = get.omittedIDs(obj)
            if isempty(obj.unreliableIDs) && isempty(obj.badSegIDs)
                out = [];
            else
                out = sort(unique([obj.unreliableIDs, obj.badSegIDs]));
            end
        end

        function classify(obj)
            % Looks at which ROIs responded to each stimulus
            obj.calcStimResponsivity();

            % Groups ROIs by which cones had significant responses
            obj.calcConeCombinations();

            % Count up cells by attribute (responsive, reliable, etc)
            obj.calcCounts();

            % Classifies different forms of opponency/non-opponency
            obj.calcOpponency();

            % Further subdivides cone groupings by cone polarity
            obj.classifyPolarity();
            % Mark cells that were "bad" and why (no response, omitted)
            obj.classifyBad();

            % Report results to the command line
            obj.reportStimResponsivity();
            obj.reportStimReliability();
            obj.reportResponseClasses();
            obj.reportOpponency();
            obj.reportPolarity();
        end

        function addOmittedIDs(obj, IDs)
            obj.omittedIDs = sort(unique(cat(2, obj.omittedIDs, IDs)));
        end

        function xy = getLMCoords(obj, coneClass)
            assert(isKey(obj.coneMap, coneClass),...
                sprintf('Key %s not recognized', coneClass));
            idx = obj.coneMap(coneClass);
            if isempty(idx)
                xy = [];
                return
            end
            x = obj.signedSNR(idx, 2);
            y = obj.signedSNR(idx, 3);
            xy = [x, y];
        end

        function out = getPctResponsive(obj, N)
            numBad = nnz(obj.isBad);
            out = round(100 * (N / (obj.numRois - numBad)), 2);
        end
    end

    methods
        function calcBadCells(obj)
            % obj.isUnresponsive = ~obj.passTable.L & ~obj.passTable.M & ~obj.passTable.S & ~obj.passTable.A;
            obj.isUnresponsive = sum(obj.passTable{:, obj.STIMULI ~= "C"}) == 0;
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

        function calcConeCombinations(obj)
            % Response classes (bad cells excluded)
            obj.responseMap = containers.Map();
            obj.responseMap('S')   =  find(~obj.passTable.L & ~obj.passTable.M &  obj.passTable.S); 
            obj.responseMap('M')   =  find(~obj.passTable.L &  obj.passTable.M & ~obj.passTable.S); 
            obj.responseMap('L')   =  find( obj.passTable.L & ~obj.passTable.M & ~obj.passTable.S); 
            obj.responseMap('A')   =  find(~obj.passTable.L & ~obj.passTable.M & ~obj.passTable.S & obj.passTable.A); 
            
            obj.responseMap('LM')  =  find( obj.passTable.L &  obj.passTable.M & ~obj.passTable.S); 
            obj.responseMap('LMS') =  find( obj.passTable.L &  obj.passTable.M &  obj.passTable.S); 
            obj.responseMap('LS')  =  find( obj.passTable.L & ~obj.passTable.M &  obj.passTable.S); 
            obj.responseMap('MS')  =  find(~obj.passTable.L &  obj.passTable.M &  obj.passTable.S); 

            obj.responseMap('NoResp')  = find(obj.isUnresponsive);
            obj.responseMap('NoResp') = setxor(obj.responseMap('NoResp'), obj.badSegIDs);
            obj.responseMap('Unrely') = obj.unreliableIDs;
            obj.responseMap('BadSeg') = obj.badSegIDs;
        end

        function calcCounts(obj)
            % Analyzed ROIs that responded significantly 
            obj.numResponsive = obj.numAnalyzed - numel(obj.responseMap('NoResp'));
            % Analyzed ROIs that responded significantly AND reliably 
            obj.numReliable = obj.numAnalyzed - numel(unique([obj.responseMap('NoResp'); obj.unreliableIDs']));
        end

        function calcOpponency(obj)
            % Helper classifiers (including bad cells)
            obj.comparisonMap = containers.Map();
            obj.comparisonMap('LvM') = find(obj.avgSign(:,2) ~= obj.avgSign(:,3));
            obj.comparisonMap('LvS') = find(obj.avgSign(:,2) ~= obj.avgSign(:,4));
            obj.comparisonMap('MvS') = find(obj.avgSign(:,3) ~= obj.avgSign(:,4));
            obj.comparisonMap('LpM') = find(obj.avgSign(:,2) == obj.avgSign(:,3)); 
            obj.comparisonMap('LpS') = find(obj.avgSign(:,2) == obj.avgSign(:,4));
            obj.comparisonMap('MpS') = find(obj.avgSign(:,3) == obj.avgSign(:,4));
        
            % Find intersections with exclusive response maps
            obj.coneMap = containers.Map();
            obj.coneMap('LvM') = intersect(obj.responseMap('LM'), obj.comparisonMap('LvM'));
            obj.coneMap('LvS') = intersect(obj.responseMap('LS'), obj.comparisonMap('LvS'));
            obj.coneMap('MvS') = intersect(obj.responseMap('MS'), obj.comparisonMap('MvS'));
            obj.coneMap('LvMS') = intersect(intersect(...
                obj.responseMap('LMS'), obj.comparisonMap('LvM')), obj.comparisonMap('MpS'));
            obj.coneMap('MvLS') = intersect(intersect(...
                obj.responseMap('LMS'), obj.comparisonMap('LvM')), obj.comparisonMap('LpS'));
            obj.coneMap('SvLM') = intersect(intersect(...
                obj.responseMap('LMS'), obj.comparisonMap('LvM')), obj.comparisonMap('LpM'));

            obj.coneMap('LS') = setdiff(obj.responseMap('LS'), obj.coneMap('LvS'));
            obj.coneMap('MS') = setdiff(obj.responseMap('MS'), obj.coneMap('MvS'));
            obj.coneMap('LM') = setdiff(obj.responseMap('LM'), obj.coneMap('LvM'));
            obj.coneMap('LMS') = intersect(intersect(...
                obj.responseMap('LMS'), obj.comparisonMap('LpM')), obj.comparisonMap('LpS'));

            obj.coneMap('S') = obj.responseMap('S');
            obj.coneMap('L') = obj.responseMap('L');
            obj.coneMap('M') = obj.responseMap('M');
        end
    end

    % Reporting methods
    methods 
        function reportStimResponsivity(obj)
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

        function reportStimReliability(obj)
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
            fprintf('L-only = %u\n', nnz(obj.responseMap('L')));
            fprintf('M-only = %u\n', nnz(obj.responseMap('M')));
            fprintf('S-only = %u\n', nnz(obj.responseMap('S')));
            fprintf('A-only = %u\n', nnz(obj.responseMap('A')));
            fprintf('LS-only = %u\n', nnz(obj.responseMap('LS')));
            fprintf('MS-only = %u\n', nnz(obj.responseMap('MS')));
            fprintf('LM-only = %u\n', nnz(obj.responseMap('LM')));
            fprintf('LMS-only = %u\n', nnz(obj.responseMap('LMS')));
            fprintf('Nonresponsive = %u\n', nnz(obj.responseMap('NoResp')));
            fprintf('Unreliable = %u\n', nnz(obj.responseMap('Unrely')));
        end

        function printConeResponse(obj)
            k = obj.responseMap.keys;
            idx = cellfun(@(x) ismember(string(x), ["NoResp", "Unrely", "BadSeg"]), k);

            fprintf('\nRESPONSE CLASSIFICATION\n');
            for i = 1:numel(k)
                if ~idx(i)
                    fprintf('\t%s = %u\n', k{i}, numel(obj.responseMap(k{i})));
                end
            end
            fprintf('\tNonresponsive = %u\n', numel(obj.responseMap('NoResp')));
            fprintf('\tUnreliable = %u\n', numel(obj.responseMap('Unrely')));
            fprintf('\tBad ROI = %u\n', nnz(obj.responseMap('BadSeg')));
        end

        function reportPolarity(obj)
            fprintf('\nPOLARITY\n');
            fprintf('L-ON = %u\n', nnz(obj.roiClasses == "L-ON"));
            fprintf('L-OFF = %u\n', nnz(obj.roiClasses == "L-OFF"));
            fprintf('M-ON = %u\n', nnz(obj.roiClasses == "M-ON"));
            fprintf('M-OFF = %u\n', nnz(obj.roiClasses == "M-OFF"));
            fprintf('S-ON = %u\n', nnz(obj.roiClasses == "S-ON"));
            fprintf('S-OFF = %u\n', nnz(obj.roiClasses == "S-OFF"));
            fprintf('A-ON = %u\n', nnz(obj.roiClasses == "A-ON"));
            fprintf('A-OFF = %u\n', nnz(obj.roiClasses == "A-OFF"));
        end
    end

    % Classification methods
    methods 
        function classifyBad(obj)
            % Cells that register as significant but only because one large
            % response is skewing the average SNR for a stimulus
            obj.roiClasses(obj.unreliableIDs) = "NR";
            % These cells did not respond significantly to any stimulus
            obj.roiClasses(obj.responseMap('NoResp')) = "NR";
            % These cells were badly segmented and should be completely 
            % omitted from the analysis
            obj.roiClasses(obj.badSegIDs) = "OMIT";
        end

        function classifyPolarity(obj)
            % Separate non-opponent cells into ON and OFF
            aOnly = obj.responseMap('A');
            for i = 1:numel(aOnly)
                if obj.avgSign(aOnly(i), end) > 0
                    obj.roiClasses(aOnly(i)) = "A-ON";
                else
                    obj.roiClasses(aOnly(i)) = "A-OFF";
                end
            end
            
            lOnly = obj.responseMap('L');
            for i = 1:numel(lOnly)
                if obj.avgSign(lOnly(i), 2) > 0
                    obj.roiClasses(lOnly(i)) = "L-ON";
                else
                    obj.roiClasses(lOnly(i)) = "L-OFF";
                end
            end
        
            mOnly = obj.responseMap('M');
            for i = 1:numel(mOnly)
                if obj.avgSign(mOnly(i), 3) > 0
                    obj.roiClasses(mOnly(i)) = "M-ON";
                else
                    obj.roiClasses(mOnly(i)) = "M-OFF";
                end
            end

            sOnly = obj.responseMap('S');
            for i = 1:numel(sOnly)
                if obj.avgSign(sOnly(i), 4) > 0
                    obj.roiClasses(sOnly(i)) = "S-ON";
                else
                    obj.roiClasses(sOnly(i)) = "S-OFF";
                end
            end

            obj.roiClasses(obj.coneMap('LvM')) = "LvM";
            obj.roiClasses(obj.coneMap('LvMS')) = "LvMS";
            obj.roiClasses(obj.coneMap('MvLS')) = "MvLS";
            obj.roiClasses(obj.coneMap('SvLM')) = "SvLM";
            
            MS = obj.coneMap('MS');
            if ~isempty(MS)
                for i = 1:numel(MS)
                    if obj.avgSign(MS(i), 4) > 0
                        obj.roiClasses(MS(i)) = "MS-ON";
                    else
                        obj.roiClasses(MS(i)) = "MS-OFF";
                    end
                end
            end

            LS = obj.coneMap('LS');
            if ~isempty(LS)
                for i = 1:numel(LS)
                    if obj.avgSign(LS(i), 4) > 0
                        obj.roiClasses(LS(i)) = "LS-ON";
                    else
                        obj.roiClasses(LS(i)) = "LS-OFF";
                    end
                end
            end

            LM = obj.coneMap('LM');
            if ~isempty(LM)
                for i = 1:numel(LM)
                    if obj.avgSign(LM(i), 2) > 0
                        obj.roiClasses(LM(i)) = "LM-ON";
                    else
                        obj.roiClasses(LM(i)) = "LM-OFF";
                    end
                end
            end
        end
    end

    %methods
    %    function plotStimResponses(obj, whichStim)
    %        stimIdx = strcmpi(obj.STIMULI, whichStim);
    %    end
    %end
end 