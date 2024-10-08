classdef Node < handle
%% The Node Class (Node)
%   Authors
%     - Bashar Tahir, btahir@nt.tuwien.ac.at
%   (c) 2016 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 

%% The Node Properties
% Accessibility
    properties
        Name            % Node name (optional)
        ID              % The unique identifier of this Node
        Type            % Node type: {'BaseStation', 'User', 'Relay'}
        PlotResults     % Indicate whether this node will get its results plotted or not.
    end
% Connectivity
    properties
        TransmitBS      % Basestations transmitting to this Node (ID).
        ReceiveBS       % Basestations receiving from this Node (ID).
        TransmitUE      % Users transmitting to this Node (ID).
        ReceiveUE       % Users receiving rom this Nodee (ID).
    end
    
%% The Node Functions
    methods
        function obj = Node(name, id, type, plotResults)
            obj.Name = name;
            obj.ID = id;
            obj.Type = type;
            obj.PlotResults = plotResults;
        end
        function attachTransmitNode(obj, txNodeID, txNodeType)      
            switch txNodeType
                case 'BaseStation'
                    obj.TransmitBS = [obj.TransmitBS txNodeID];
                case 'User'
                    obj.TransmitUE = [obj.TransmitUE txNodeID];
            end
        end
        function attachReceiveNode(obj, rxNodeID, rxNodeType)      
            switch rxNodeType
                case 'BaseStation'
                    obj.ReceiveBS = [obj.ReceiveBS rxNodeID];
                case 'User'
                    obj.ReceiveUE = [obj.ReceiveUE rxNodeID];
            end
        end
     end   
end

