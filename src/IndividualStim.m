classdef IndividualStim < ConeDataset

    properties 
        STIMULI = ["C", "L", "M", "S", "A"];
    end

    methods
        function obj = IndividualStim(animalID, dataset, varargin)
            obj@ConeDataset(animalID, dataset, varargin{:});
        end
    end

    % Subclass-specific analyses
    methods
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
    end

    methods (Access = protected)
        function calcConeCombinations(obj)
            calcConeCombinations@ConeDataset(obj);
            
            
            % Response classes (bad cells excluded)
            obj.responseMap('S')   =  find(~obj.passTable.L & ~obj.passTable.M &  obj.passTable.S); 
            obj.responseMap('M')   =  find(~obj.passTable.L &  obj.passTable.M & ~obj.passTable.S); 
            obj.responseMap('L')   =  find( obj.passTable.L & ~obj.passTable.M & ~obj.passTable.S); 
            obj.responseMap('A')   =  find(~obj.passTable.L & ~obj.passTable.M & ~obj.passTable.S & obj.passTable.A); 
            
            obj.responseMap('LM')  =  find( obj.passTable.L &  obj.passTable.M & ~obj.passTable.S); 
            obj.responseMap('LMS') =  find( obj.passTable.L &  obj.passTable.M &  obj.passTable.S); 
            obj.responseMap('LS')  =  find( obj.passTable.L & ~obj.passTable.M &  obj.passTable.S); 
            obj.responseMap('MS')  =  find(~obj.passTable.L &  obj.passTable.M &  obj.passTable.S); 
        end
        
        function calcComparisons(obj)
            % Classify cone comparisons (not exclusive)
            obj.comparisonMap('LvM') = find(obj.avgSign(:,2) ~= obj.avgSign(:,3));
            obj.comparisonMap('LvS') = find(obj.avgSign(:,2) ~= obj.avgSign(:,4));
            obj.comparisonMap('MvS') = find(obj.avgSign(:,3) ~= obj.avgSign(:,4));
            obj.comparisonMap('LpM') = find(obj.avgSign(:,2) == obj.avgSign(:,3)); 
            obj.comparisonMap('LpS') = find(obj.avgSign(:,2) == obj.avgSign(:,4));
            obj.comparisonMap('MpS') = find(obj.avgSign(:,3) == obj.avgSign(:,4));
        end

        function calcOpponency(obj)
            % Single cone input cells
            obj.coneMap('S') = obj.responseMap('S');
            obj.coneMap('L') = obj.responseMap('L');
            obj.coneMap('M') = obj.responseMap('M');

            % Multiple cone input cells
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
        end
    end 
end 