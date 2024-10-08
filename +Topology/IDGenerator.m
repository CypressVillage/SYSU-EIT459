classdef IDGenerator < handle
%% The IDGenerator Class (IDGenerator)
%   Authors
%     - Bashar Tahir, btahir@nt.tuwien.ac.at
%   (c) 2016 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 

%% The IDGenerator Properties
    properties
        currentID = 0;
    end 
%% The IDGenerator Functions
    methods
        function obj = IDGenerator()
        end
        function newID = generateID(obj)
            % At the moment, it is just a simple addition-by-one rule.
            obj.currentID = obj.currentID + 1;
            newID = obj.currentID;
        end
    end   
    
end

