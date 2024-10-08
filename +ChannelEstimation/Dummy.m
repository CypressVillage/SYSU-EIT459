classdef Dummy < handle
    % dummy channel estimation object for perfect channel estimation
    
    properties
        PilotMatrix
    end
    
    methods
        function obj = Dummy( nSubcarriers, nSymbols )
            % dummy class constructor
            
            obj.PilotMatrix = zeros( nSubcarriers, nSymbols );
        end
    end
    
end

