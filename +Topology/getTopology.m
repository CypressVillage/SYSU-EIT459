%% Topology Generation
%   Authors
%     - Bashar Tahir, btahir@nt.tuwien.ac.at
%   (c) 2018 Institute of Telecommunications, TU Wien.
%   www.nt.tuwien.ac.at 
%%
function [Links, BS, UE] = getTopology(simParams)
% Parse simParams.topology, and construct the corresponding network/links
    BSnum = 0;
    UEnum = 0;
    %/ Checks
    if length(simParams.topology.nodes) < 7
        error('Alteast two nodes have to be specified in the topology.');
    end
    %/
    IDGenerator = Topology.IDGenerator;
    split1 = strsplit(simParams.topology.nodes(1:end),',');
    for ii = 1:length(split1)
        index = str2double(split1{ii}(3:end));
        %/ Checks
         if isnan(index)
             error('The topology nodes are entered incorrectly. Please make sure no extra characters or spaces are present.')
         end
          if index == 0
             error('Node index of zero is not allowed (i.e. BS0 or UE0).')
         end
        %/
        switch split1{ii}(1:2)
            case 'BS'
                if index == BSnum + 1
                    ID = IDGenerator.generateID();
                    if length(simParams.simulation.plotResultsFor) == 1
                        plotResults = simParams.simulation.plotResultsFor;
                    else
                        plotResults = simParams.simulation.plotResultsFor(ID);
                    end
                    BS{index} = Elements.BaseStation(['BS' num2str(index)], ID, simParams.simulation.nAntennasBaseStation(index), simParams.simulation.txPowerBaseStation(index), plotResults); %#ok
                    BSnum = index;
                %/ Checks
                elseif index == BSnum
                    error(['BS' num2str(index) ' already exists. Please check topology.nodes for duplicates.'])
                else
                    error('The nodes indices are not entered in an ascending order. Make sure to have them in the correct order, for example: BS1,BS2,BS3,UE1,UE2,UE3.')
                end
                %/
            case 'UE'
                if index == UEnum + 1
                    ID = IDGenerator.generateID();
                    if length(simParams.simulation.plotResultsFor) == 1
                        plotResults = simParams.simulation.plotResultsFor;
                    else
                        plotResults = simParams.simulation.plotResultsFor(ID);
                    end
                    UE{index} = Elements.User(['UE' num2str(index)], ID, simParams.simulation.nAntennasUser(index), simParams.simulation.txPowerUser(index), plotResults); %#ok
                    UEnum = index;
                %/ Checks
                elseif index == UEnum
                    error(['UE' num2str(index) ' already exists. Please check topology.nodes for duplicates.'])
                else
                    error('The nodes indices are not entered in an ascending order. Make sure to have them in the correct order, for example: BS1,BS2,BS3,UE1,UE2,UE3.')    
                end
                %/
            otherwise
                error('The nodes name must start with either "BS" or "UE".');
        end
    end
    Links = cell(length(split1));
    
    % Generation of Primary Links
    split1 = strsplit(simParams.topology.primaryLinks,',');
    linkID = 1;
    nUplinks = 0;
    nDownlinks = 0;
    nD2Dlinks = 0;
    for ii = 1:length(split1)
        split2 = strsplit(split1{ii},':');
        if strcmp(split2, '') == 0
            % transmitter index
            txIndex = str2double(split2{1}(3:end));
            
            %/ Checks
            if isnan(txIndex)
                error('The nodes in topology.primaryLinks are entered incorrectly. Please make sure no extra characters or spaces are present.')
            end
            %/
            switch split2{1}(1:2)
                case 'BS'
                    %/ Checks
                    if txIndex > BSnum
                        error(['BS' num2str(txIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                    end
                    %/
                    
                    tx              = BS{txIndex}; 
                case 'UE'
                    %/ Checks
                    if txIndex > UEnum
                        error(['UE' num2str(txIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                    end
                    %/
                    tx              = UE{txIndex};
                    ueIndex         = txIndex;
                otherwise
                    error('The nodes name must start with either "BS" or "UE". Please check the entered primaryLinks.');
            end
            
            % receiver index
            rxIndex = str2double(split2{2}(3:end));
            %/ Checks
            if isnan(rxIndex)
                error('The nodes in topology.primaryLinks are entered incorrectly. Please make sure no extra characters or spaces are present.')
            end
            %/
            switch split2{2}(1:2)
                case 'BS'
                    %/ Checks
                    if rxIndex > BSnum
                        error(['BS' num2str(rxIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                    end
                    %/
                    
                    rx      = BS{rxIndex};
                case 'UE'
                    %/ Checks
                    if rxIndex > UEnum
                        error(['UE' num2str(rxIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                    end
                    %/
                    
                    rx      = UE{rxIndex};
                    ueIndex = rxIndex;
                otherwise
                    error('The nodes name must start with either "BS" or "UE". Please check the entered primaryLinks.');
            end              
            rx.attachTransmitNode(tx.ID, tx.Type);
            tx.attachReceiveNode(rx.ID, rx.Type);
            
            % assign direction for the link
            if(strcmp(tx.Type,'BaseStation') && strcmp(rx.Type,'User'))
                direction = 'Downlink';
                nDownlinks = nDownlinks + 1;
            elseif(strcmp(tx.Type,'User') && strcmp(rx.Type,'BaseStation'))
                direction = 'Uplink';
                nUplinks = nUplinks + 1;
            elseif(strcmp(tx.Type,'User') && strcmp(rx.Type,'User'))
                direction = 'D2D';
                nD2Dlinks = nD2Dlinks + 1;
            else
                direction = [];
            end
            
            %/ Checks
            if ~isempty(Links{tx.ID, rx.ID})
                error(['A link from the ' tx.Type ' with ID ' num2str(tx.ID)  ' to the ' rx.Type ' with ID ' num2str(rx.ID) ', already exists. Please check for any duplicate links in topology.primaryLinks.']);
            end
            %/

            Links{tx.ID, rx.ID} =  Elements.Link(linkID, 'Primary', direction, tx.ID, rx.ID, 'CC', simParams.simulation.userVelocity(ueIndex) , 0, tx.nAntennas, rx.nAntennas,simParams.feedback,simParams.modulation,simParams.modulation.transmissionMode);        
            linkID = linkID + 1;
        end
    end
    %/ Checks
    if nDownlinks == 0 && simParams.simulation.simulateDownlink
        error('Downlink simulations are selected, however, there is not a single link from a user to a base station. To simulate downlink, add at least one downlink to the topology.');
    end
    if nUplinks == 0 && simParams.simulation.simulateUplink
        error('Uplink simulations are selected, however, there is not a single link from a base station to a user. To simulate uplink, add at least one uplink to the topology.');
    end
    if nD2Dlinks == 0 && simParams.simulation.simulateD2D
        error('D2D simulations are selected, however, there is not a single link from a user to a user. To simulate D2D, add at least one D2D link to the topology.');
    end
    %/
    
    % Generation of Interference Links
    switch simParams.topology.interferenceGeneration
        case 'Automatic'
            %/ Checks
            if simParams.topology.attenuation < 0
                error('The attenuation value of the interference links has to be positive. Please check topology.attenuation.');
            end
            %/
            for iBS = 1:length(BS)
                BSID = BS{iBS}.ID;
                for iUE = 1:length(UE)
                    UEID = UE{iUE}.ID;
                    if isempty(Links{BSID, UEID}) && simParams.simulation.simulateDownlink
                        UE{iUE}.attachTransmitNode(BSID, 'BaseStation');
                        BS{iBS}.attachReceiveNode(UEID, 'User');

                        Links{BSID, UEID} =  Elements.Link(0, 'Interference', 'Downlink', BSID, UEID, 'CC', simParams.simulation.userVelocity(iUE), simParams.topology.attenuation, BS{iBS}.nAntennas, UE{iUE}.nAntennas,simParams.feedback,simParams.modulation,simParams.modulation.transmissionMode);

                    end   
                    if isempty(Links{UEID, BSID}) && simParams.simulation.simulateUplink
                        BS{iBS}.attachTransmitNode(UEID, 'User');
                        UE{iUE}.attachReceiveNode(BSID, 'BaseStation');

                        Links{UEID, BSID} =  Elements.Link(0, 'Interference', 'Uplink', UEID, BSID, 'CC', simParams.simulation.userVelocity(iUE), simParams.topology.attenuation, UE{iUE}.nAntennas, BS{iBS}.nAntennas,simParams.feedback,simParams.modulation,simParams.modulation.transmissionMode);

                    end  
                end
            end
        case 'Manual'
            split1 = strsplit(simParams.topology.interferingLinks,',');
            for ii = 1:length(split1)
                split2 = strsplit(split1{ii},':');
                if strcmp(split2, '') == 0
                    % transmitter index
                    txIndex = str2double(split2{1}(3:end));
                    %\ Checks
                    if isnan(txIndex)
                        error('The nodes in topology.interferingLinks are entered incorrectly. Please make sure no extra characters or spaces are present.')
                    end
                    %\
                    switch split2{1}(1:2)
                        case 'BS'
                            %/ Checks
                            if txIndex > BSnum
                                error(['BS' num2str(txIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                            end
                            %/
                            tx              = BS{txIndex};
                        case 'UE'
                            %/ Checks
                            if txIndex > UEnum
                                error(['UE' num2str(txIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                            end
                            %/
                            tx              = UE{txIndex};
                            ueIndex         = rxIndex;
                    end
                    
                    % receiver index
                    split3 = strsplit(split2{2},'*');
                    rxIndex = str2double(split3{1}(3:end));
                    
                    %\ Checks
                    if isnan(rxIndex)
                        error('The nodes in topology.interferingLinks are entered incorrectly. Please make sure no extra characters or spaces are present.')
                    end
                    %\
                    switch split3{1}(1:2)
                        case 'BS'
                            %/ Checks
                            if rxIndex > BSnum
                                error(['BS' num2str(rxIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                            end
                            %/
                            
                            rx      = BS{rxIndex};   
                        case 'UE'
                            %/ Checks
                            if rxIndex > UEnum
                                error(['UE' num2str(rxIndex) ' does not exist. Make sure it is added in topology.nodes.']); 
                            end
                            %/
                            
                            rx      = UE{rxIndex};
                            ueIndex = rxIndex;
                    end              
                    % assign direction for the link
                    if(strcmp(tx.Type,'BaseStation') && strcmp(rx.Type,'User'))
                        direction = 'Downlink';
                    elseif(strcmp(tx.Type,'User') && strcmp(rx.Type,'BaseStation'))
                        direction = 'Uplink';
                    elseif(strcmp(tx.Type,'User') && strcmp(rx.Type,'User'))
                        direction = 'D2D';
                    else
                        direction = [];
                    end
                    rx.attachTransmitNode(tx.ID, tx.Type);
                    tx.attachReceiveNode(rx.ID, rx.Type);
                    
                    %/ Checks
                    if ~isempty(Links{tx.ID, rx.ID})
                        error(['A link from the ' tx.Type ' with ID ' num2str(tx.ID)  ' to the ' rx.Type ' with ID ' num2str(rx.ID) ', already exists. Please check for any duplicate links in topology.interferingLinks.']);
                    end             
                    if str2double(split3{2}) < 0
                        error('The attenuation value of the interference links has to be positive.');
                    end
                    %/                

                    Links{tx.ID, rx.ID} =  Elements.Link(0, 'Interference', direction, tx.ID, rx.ID, 'CC', simParams.simulation.userVelocity(ueIndex), str2double(split3{2}), tx.nAntennas, rx.nAntennas,simParams.feedback,simParams.modulation,simParams.modulation.transmissionMode);       

                end
            end
        otherwise
            error('Interference generation type is invalid. Please check topology.interferenceGeneration.');
    end
    
    %/ Checks
      for iBS = 1:length(BS)
        if isempty(BS{iBS}.TransmitBS) && isempty(BS{iBS}.ReceiveBS) && isempty(BS{iBS}.TransmitUE) && isempty(BS{iBS}.ReceiveUE)
            error(['The Basestation with ID ' num2str(BS{iBS}.ID) ' is not connected to any other node.'])
        end
      end
      for iUE = 1:length(UE)
        if isempty(UE{iUE}.TransmitBS) && isempty(UE{iUE}.ReceiveBS) && isempty(UE{iUE}.TransmitUE) && isempty(UE{iUE}.ReceiveUE)
            error(['The user with ID ' num2str(UE{iUE}.ID) ' is not connected to any other node.'])
        end
      end
    %/
end