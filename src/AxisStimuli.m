classdef AxisStimuli < ConeDataset

    properties (Hidden)
        STIMULI = ["I", "S", "A"];
        HUETUNED = ["L-MS", "M-LS", "MS-L", "LS-M"];
    end

    methods 
        function obj = AxisStimuli(animalID, dataset, varargin)
            obj@ConeDataset(animalID, dataset, varargin{:});
        end
    end

    methods (Access = protected)
        function calcConeCombinations(obj)
            calcConeCombinations@ConeDataset(obj);
            % Response classes (bad cells excluded)
            obj.responseMap('S')   =  find(~obj.passTable.I &  obj.passTable.S); 
            obj.responseMap('I')   =  find( obj.passTable.I & ~obj.passTable.S); 
            obj.responseMap('A')   =  find(~obj.passTable.I & ~obj.passTable.S & obj.passTable.A);

            obj.responseMap('SI')  =  find( obj.passTable.I &  obj.passTable.S);
        end

        function calcComparisons(obj)
            obj.comparisonMap('IvS') = intersect(obj.responseMap('SI'),... 
                find(obj.avgSign(:, 1) ~= obj.avgSign(:, 2)));
            obj.comparisonMap('IpS') = intersect(obj.responseMap('SI'), ...
                find(obj.avgSign(:, 1) == obj.avgSign(:, 2)));
            obj.comparisonMap('LvM') = intersect(obj.responseMap('I'), find(obj.avgSign(:,1) > 0));
            obj.comparisonMap('MvL') = setdiff(obj.responseMap('I'), obj.comparisonMap('LvM'));
            obj.comparisonMap('LvMaS') = intersect(obj.responseMap('SI'), find(obj.avgSign(:,1) > 0));
            obj.comparisonMap('MvLaS') = intersect(obj.responseMap('SI'), find(obj.avgSign(:,1) < 0));
        end

        function calcOpponency(obj)
            obj.coneMap('S') = obj.responseMap('S');
            obj.coneMap('I') = obj.responseMap('I');
            obj.coneMap('LvM') = setdiff(obj.comparisonMap('LvM'), obj.responseMap('S'));
            obj.coneMap('MvL') = setdiff(obj.comparisonMap('MvL'), obj.responseMap('S'));
            % L-M (+), LvMvS (-) = L-M-S = L-MS
            obj.coneMap('L-MS') = intersect(obj.comparisonMap('LvMaS'), obj.comparisonMap('IvS'));
            % M-L (-), LvMvS (+) = 
            obj.coneMap('MS-L') = intersect(obj.comparisonMap('MvLaS'), obj.comparisonMap('IvS'));
            % L-M (+), LvMpS (+) = L-M+S = LS-M
            obj.coneMap('LS-M') = intersect(obj.comparisonMap('LvMaS'), obj.comparisonMap('IpS'));
            % M-L (-), LvMpS (-) = M-L-S = M-LS
            obj.coneMap('M-LS') = intersect(obj.comparisonMap('MvLaS'), obj.comparisonMap('IpS'))
        end
    end

    % Extended methods
    methods 
        function classifyPolarity(obj)
            classifyPolarity@ConeDataset(obj);

            % Remove the ON/OFF stuff for L-M and S comparisons
            for i = 1:numel(obj.HUETUNED)
                iRoi = obj.coneMap(char(obj.HUETUNED(i)));
                if ~isempty(iRoi)
                    obj.roiClasses(iRoi) = obj.HUETUNED(i);
                end
            end
        end
    end
end