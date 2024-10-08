function [ codebook] = getCodebook( nTx , transmissionMode, varargin)

switch transmissionMode
    
    case 'CLSM'
        if nTx~=1 &&  nTx~=2 && nTx~=4 && nTx~=8
            error('The transmission mode %s only supports 2, 4 or 8 transmit antennas',transmissionMode);
        end
            switch nTx
                case 1
                    codebook.W{1}=1;
                case 2
                    codebook.W = cell(7);
                    codebook.W{1}=[1;1]/sqrt(2);
                    codebook.W{2}=[1;-1]/sqrt(2);
                    codebook.W{3}=[1;1i]/sqrt(2);
                    codebook.W{4}=[1;-1j]/sqrt(2);
                    codebook.W{5}=[1,0;0,1]/sqrt(2);
                    codebook.W{6}=[1,1;1,-1]/2;
                    codebook.W{7}=[1 1;1j,-1j]/2;
                case 4

                    V = [  1,-1,-1,-1;     % Matrix corresponding to vectors u0 ... u15 in Table 6.3.4.2.3-2
                        1,-1i,1,1i;
                        1,1,-1,1;
                        1,1i,1,-1i;
                        1,(-1-1i)/sqrt(2), -1i,(1-1i)/sqrt(2);
                        1,(1-1i)/sqrt(2), 1i,(-1-1i)/sqrt(2);
                        1,(1+1i)/sqrt(2), -1i,(-1+1i)/sqrt(2);
                        1,(-1+1i)/sqrt(2), 1i,(1+1i)/sqrt(2);
                        1,-1,1,1;
                        1,-1i,-1,-1i;
                        1,1,1,-1;
                        1,1i,-1,1i;
                        1,-1,-1,1;
                        1,-1,1,-1;
                        1,1,-1,-1;
                        1,1,1,1;].';

                    W=zeros(4,4,16);
                    I=eye(4,4);

                    for i=1:16
                        W(:,:,i)=I-2*V(:,i)*V(:,i)'/(V(:,i)'*V(:,i));
                    end

                    LTEmapping{1} = ones(16,1);
                    LTEmapping{2} = [1 4; 1 2; 1 2; 1 2; 1 4; 1 4; 1 3; 1 3; 1 2; 1 4; 1 3; 1 3; 1 2; 1 3; 1 3; 1 2];
                    LTEmapping{3} = [1 2 4; 1 2 3; 1 2 3; 1 2 3; 1 2 4; 1 2 4; 1 3 4; 1 3 4; 1 2 4; 1 3 4; 1 2 3; 1 3 4; 1 2 3; 1 2 3; 1 2 3; 1 2 3];
                    LTEmapping{4} = [1 2 3 4; 1 2 3 4; 3 2 1 4; 3 2 1 4; 1 2 3 4; 1 2 3 4; 1 3 2 4; 1 3 2 4;
                        1 2 3 4; 1 2 3 4; 1 3 2 4; 1 3 2 4; 1 2 3 4; 1 3 2 4; 3 2 1 4; 1 2 3 4];
                   codebook.W=cell(1,16*4);
                    for nLayers=1:4
                        temp_W = 1/sqrt(nLayers)*W;
                        for codebook_index=0:15
                            codebook.W{(nLayers-1)*16+codebook_index+1} = temp_W(:,LTEmapping{nLayers}(codebook_index+1,:),codebook_index+1);
                        end
                    end
                case 8
                    codebook.W = codebook_8TX;
                otherwise
                    error('No codebook for %d transmit antennas',nTx);
            end
            codebook.U = 0;
            codebook.D = 0;
            codebook.methode = 'CLSM';
    case 'OLSM'
        if nTx~=1 &&  nTx~=2 && nTx~=4
            error('The transmission mode %s only supports 2 or 4 transmit antennas',transmissionMode);
        end
        switch nTx
            case 1
                codebook.W{1}=1;
                 
            case 2   % precoding matrix according to Table 6.3.4.2.3-1
                codebook.W = cell(1,2);
                codebook_index = 1;
                switch codebook_index
                    case 0
                        
                           codebook.W{1,2} = 1/sqrt(2)*[1,0;0,1];
                           codebook.W{1,1} = 1/sqrt(2)*[1;1];
                        
                    case 1
                        
                            codebook.W{1,2} = 1/(2)*[1,1;1,-1];
                            codebook.W{1,1} = 1/sqrt(2)*[1;-1];
                       
                    case 2
                       
                            codebook.W{1,2} = 1/(2)*[1,1;1i,-1i];
                            codebook.W{1,1} = 1/sqrt(2)*[1;1i];
                       
                    case 3
                        
                            codebook.W{1,2} = 0;
                            codebook.W{1,1} = 1/sqrt(2)*[1;-1];
                        
                    otherwise
                        error('codebook index not supported');
                end
            case 4
                
                codebook_index = [12,13,14,15];
                codebook.W = cell(size(codebook_index,2),nTx);
                for i = 1:length(codebook_index)    % for open loop spatial multiplexing 
                                                        % codebook_index = [12,13,14,15]
                                                        % see 3GPP TS 36.213-820, Section 7.1.3 
                                                        % NOTE: store this in some pre-stored variable
                        switch codebook_index(i)
                            case 0
                                u= [1;-1;-1;-1];
                            case 1
                                u= [1;-1i;1;1i];
                            case 2
                                u= [1;1;-1;1];
                            case 3
                                u= [1;1i;1;-1i];
                            case 4
                                u= [1;(-1-1i)/sqrt(2); -1i;(1-1i)/sqrt(2)];
                            case 5
                                u= [1;(1-1i)/sqrt(2); 1i;(-1-1i)/sqrt(2)];
                            case 6
                                u= [1;(1+1i)/sqrt(2); -1i;(-1+1i)/sqrt(2)];
                            case 7
                                u= [1;(-1+1i)/sqrt(2); 1i;(1+1i)/sqrt(2)];
                            case 8
                                u= [1;-1;1;1];
                            case 9
                                u= [1;-1i;-1;-1i];
                            case 10
                                u= [1;1;1;-1];
                            case 11
                                u= [1;1i;-1;1i];
                            case 12
                                u= [1;-1;-1;1];
                            case 13
                                u= [1;-1;1;-1];
                            case 14
                                u= [1;1;-1;-1];
                            case 15
                                u= [1;1;1;1];
                            otherwise
                                error('codebook index not supported');
                        end
                        W1 = diag(ones(1,4))-2*u*u'/(u'*u);
                        for RI = 1:nTx
                                switch RI  
                                case 1
                                    codebook.W{i,RI} = W1(:,1);
                                case 2
                                    switch codebook_index(i)
                                        case {0,4,5,9}
                                            codebook.W{i,RI} = W1(:,[1 4])/sqrt(2);
                                        case {1,2,3,8,12,15}
                                            codebook.W{i,RI} = W1(:,[1 2])/sqrt(2);
                                        otherwise
                                            codebook.W{i,RI} = W1(:,[1 3])/sqrt(2);
                                    end
                                case 3
                                    switch codebook_index(i)
                                        case {0,4,5,8}
                                            codebook.W{i,RI} = W1(:,[1 2 4])/sqrt(3);
                                        case {1,2,3,10,12,13,14,15}
                                            codebook.W{i,RI} = W1(:,[1 2 3])/sqrt(3);
                                        otherwise
                                            codebook.W{i,RI} = W1(:,[1 3 4])/sqrt(3);
                                    end
                                case 4
                                    switch codebook_index(i)
                                        case {0,1,4,5,8,9,12,15}
                                            codebook.W{i,RI} = W1(:,[1 2 3 4])/2;
                                        case {6,7,10,11,13}
                                            codebook.W{i,RI} = W1(:,[1 3 2 4])/2;
                                        otherwise
                                            codebook.W{i,RI} = W1(:,[3 2 1 4])/2;
                                    end
                                otherwise
                                    error('RI not supported');
                                end
                        end
                end
            otherwise
                error('No codebook for %d transmit antennas',nTx);
        end
        
        
        if length(varargin)<1
            CDD = 0; % default mode?
        else
            CDD = varargin{1};
        end
        switch CDD  % cyclic delay diversity, precoding according to Section 6.3.4.2
            case 0  % zero delay CDD
                codebook.D{1} = 1;  
                codebook.U{1} = 1;
                codebook.D{2} = 1;  
                codebook.U{2} = 1;
                codebook.D{3} = 1;  
                codebook.U{3} = 1;
                codebook.D{4} = 1;  
                codebook.U{4} = 1;
         
            case 1  % large delay CDD
                % RI = 1
                    codebook.D{1} = 1;  
                    codebook.U{1} = 1;
                 % RI = 2
                    codebook.U{2} = 1/sqrt(2)*[1,1;1,exp(-1i*pi)];
                    codebook.D{2} = [1,0;0,exp(-1i*pi)];
                % RI = 3
                   codebook.U{3} = 1/sqrt(3)*[1,1,1;1,exp(-1i*2*pi/3),exp(-1i*4*pi/3);1,exp(-1i*4*pi/3),exp(-1i*8*pi/3)];
                    codebook.D{3} = [1,0,0;0,exp(-1i*2*pi/3),0;0,0,exp(-1i*4*pi/3)];
                % RI = 4
                    codebook.U{4} = 1/2*[1,1,1,1;1,exp(-1i*2*pi/4),exp(-1i*4*pi/4),exp(-1i*6*pi/4);...
                        1,exp(-1i*4*pi/4),exp(-1i*8*pi/4),exp(-1i*12*pi/4);...
                        1,exp(-1i*6*pi/4),exp(-1i*12*pi/4),exp(-1i*18*pi/4)];
                    codebook.D{4} = [1,0,0,0;0,exp(-1i*2*pi/4),0,0;0,0,exp(-1i*4*pi/4),0;0,0,0,exp(-1i*6*pi/4)];
                    
            otherwise
                    error('RI not supported');

        end
        codebook.methode = 'OLSM';
        
    case 'TxD'
        % Add codebook for transmit diversity
         if nTx~=1 &&  nTx~=2 && nTx~=4 && nTx~=8
            error('The transmission mode %s only supports 2 or 4 transmit antennas',transmissionMode);
        end
        codebook.Z{1} =  [1, 0, 1i,  0;
         0,-1,  0, 1i;
         0, 1,  0, 1i;
         1, 0,-1i,  0];
     
        codebook.Z{2} =  [1, 0, 0, 0, 1i,  0,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0,-1, 0, 0,  0, 1i,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 1, 0, 0,  0, 1i,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         1, 0, 0, 0,-1i,  0,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 1, 0,  0,  0, 1i, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 0,-1,  0,  0,  0,1i;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 0, 1,  0,  0,  0,1i;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 1, 0,  0,  0,-1i, 0];
     
        codebook.Z{3} = zeros(8*8,16);
        % 1st Alamouti scheme
        codebook.Z{3}(1,1) = 1;
        codebook.Z{3}(1,9) = 1i;
        codebook.Z{3}(5,2) = -1;
        codebook.Z{3}(5,10) = 1i;
        codebook.Z{3}(9,2) = 1;
        codebook.Z{3}(9,10) = 1i;
        codebook.Z{3}(13,1) = 1;
        codebook.Z{3}(13,9) = -1i;
        % 2nd Alamouti scheme
        codebook.Z{3}(1+17,1+2) = 1;
        codebook.Z{3}(1+17,9+2) = 1i;
        codebook.Z{3}(5+17,2+2) = -1;
        codebook.Z{3}(5+17,10+2) = 1i;
        codebook.Z{3}(9+17,2+2) = 1;
        codebook.Z{3}(9+17,10+2) = 1i;
        codebook.Z{3}(13+17,1+2) = 1;
        codebook.Z{3}(13+17,9+2) = -1i;
        % 3rd Alamouti scheme
        codebook.Z{3}(1+17*2,1+2*2) = 1;
        codebook.Z{3}(1+17*2,9+2*2) = 1i;
        codebook.Z{3}(5+17*2,2+2*2) = -1;
        codebook.Z{3}(5+17*2,10+2*2) = 1i;
        codebook.Z{3}(9+17*2,2+2*2) = 1;
        codebook.Z{3}(9+17*2,10+2*2) = 1i;
        codebook.Z{3}(13+17*2,1+2*2) = 1;
        codebook.Z{3}(13+17*2,9+2*2) = -1i;
        % 4th Alamouti scheme
        codebook.Z{3}(1+17*3,1+2*3) = 1;
        codebook.Z{3}(1+17*3,9+2*3) = 1i;
        codebook.Z{3}(5+17*3,2+2*3) = -1;
        codebook.Z{3}(5+17*3,10+2*3) = 1i;
        codebook.Z{3}(9+17*3,2+2*3) = 1;
        codebook.Z{3}(9+17*3,10+2*3) = 1i;
        codebook.Z{3}(13+17*3,1+2*3) = 1;
        codebook.Z{3}(13+17*3,9+2*3) = -1i;
     
        codebook.methode = 'TxD';
    otherwise
        error('Transmission mode not supported');


end
end

function W = codebook_8TX   
    W1_back = cell(8,1);
    W2_back = cell(8,1);
    %% Rank 1
    B = zeros(4,32);
    for m = 0:3
        for n = 0:31
            B(m+1,n+1) = exp(1i*2*pi*m*n/32);
        end
    end
    X = cell(16,1);
    for k = 0:15
        X{k+1} = 1/2*[B(:,mod(2*k,32)+1),B(:,mod(2*k+1,32)+1),B(:,mod(2*k+2,32)+1),B(:,mod(2*k+3,32)+1)];
    end
    W1 = zeros(8,8,16);
    for k = 0:15
        W1(:,:,k+1) = [X{k+1},zeros(4,4);zeros(4,4),X{k+1}];
    end
    Y = eye(4);
    W2 = zeros(8,1,16);
    for k1 = 1:4
        W2(:,:,(k1-1)*4+1) = 1/sqrt(2)*[Y(:,k1);Y(:,k1)];
        W2(:,:,(k1-1)*4+2) = 1/sqrt(2)*[Y(:,k1);1i*Y(:,k1)];
        W2(:,:,(k1-1)*4+3) = 1/sqrt(2)*[Y(:,k1);-Y(:,k1)];
        W2(:,:,(k1-1)*4+4) = 1/sqrt(2)*[Y(:,k1);-1i*Y(:,k1)];
    end
    W1_back{1} = W1;
    W2_back{1} = W2;

    %% Rank 2
    W2 = zeros(8,2,16);
    Y = zeros(4,16);
    Y([1],[1,2,9,13]) = 1;
    Y([2],[3,4,10,11,15]) = 1;
    Y([3],[5,6,12]) = 1;
    Y([4],[7,8,14,16]) = 1;
    for k1 = 1:8
        for k2 = 1:2
            if k2 == 1
                W2(:,:,2*(k1-1)+1) = 1/sqrt(2)*[Y(:,(k1-1)*2+1),Y(:,(k1-1)*2+2);Y(:,(k1-1)*2+1),-Y(:,(k1-1)*2+2)];
            else
                W2(:,:,2*(k1-1)+2) = 1/sqrt(2)*[Y(:,(k1-1)*2+1),Y(:,(k1-1)*2+2);1i*Y(:,(k1-1)*2+1),-1i*Y(:,(k1-1)*2+2)];
            end
        end
    end
    W1_back{2} = W1;
    W2_back{2} = W2;    
    
    %% Rank 3
    B = zeros(4,16);
    for m = 0:3
        for n = 0:15
            B(m+1,n+1) = exp(1i*2*pi*m*n/16);
        end
    end
    X = cell(4,1);
    for k = 0:3
        X{k+1} = 1/2*[B(:,mod(4*k,16)+1),B(:,mod(4*k+1,16)+1),B(:,mod(4*k+2,16)+1),B(:,mod(4*k+3,16)+1),B(:,mod(4*k+4,16)+1),B(:,mod(4*k+5,16)+1),B(:,mod(4*k+6,16)+1),B(:,mod(4*k+7,16)+1)];
    end
    W1 = zeros(8,16,4);
    for k = 0:3
        W1(:,:,k+1) = [X{k+1},zeros(4,8);zeros(4,8),X{k+1}];
    end
    W2 = zeros(16,3,16);
    for k1 = 1:16
        if k1 <= 8
            y1 = zeros(8,1);
            y2 = zeros(8,2);
        else
            y1 = zeros(8,2);
            y2 = zeros(8,1);
        end
        switch k1
            case 1
                y1(1) = 1;
                y2(1,1) = 1;
                y2(5,2) = 1;
            case 2
                y1(2) = 1;
                y2(2,1) = 1;
                y2(6,2) = 1;
            case 3
                y1(3) = 1;
                y2(3,1) = 1;
                y2(7,2) = 1;
            case 4
                y1(4) = 1;
                y2(4,1) = 1;
                y2(8,2) = 1;
            case 5
                y1(5) = 1;
                y2(1,1) = 1;
                y2(5,2) = 1;
            case 6
                y1(6) = 1;
                y2(2,1) = 1;
                y2(6,2) = 1;
            case 7
                y1(7) = 1;
                y2(3,1) = 1;
                y2(7,2) = 1;
            case 8
                y1(8) = 1;
                y2(4,1) = 1;
                y2(8,2) = 1;
            case 9
                y1(1,1) = 1;
                y1(5,2) = 1;
                y2(5) = 1;
            case 10
                y1(2,1) = 1;
                y1(6,2) = 1;
                y2(6) = 1;
            case 11
                y1(3,1) = 1;
                y1(7,2) = 1;
                y2(7) = 1;
            case 12
                y1(4,1) = 1;
                y1(8,2) = 1;
                y2(8) = 1;
            case 13
                y1(5,1) = 1;
                y1(1,2) = 1;
                y2(1) = 1;
            case 14
                y1(6,1) = 1;
                y1(2,2) = 1;
                y2(2) = 1;
            case 15
                y1(7,1) = 1;
                y1(3,2) = 1;
                y2(3) = 1;
            case 16
                y1(8,1) = 1;
                y1(4,2) = 1;
                y2(4) = 1;
        end
        W2(:,:,k1) = 1/sqrt(2)*[y1,y2;y1,-y2];
    end
    W1_back{3} = W1;
    W2_back{3} = W2;
    %% Rank 4
    W2 = zeros(16,4,8);
    for k1 = 1:4
        Y = zeros(8,2);
        switch k1
            case 1
                Y(1,1) = 1;
                Y(5,2) = 1;
            case 2
                Y(2,1) = 1;
                Y(6,2) = 1;
            case 3
                Y(3,1) = 1;
                Y(7,2) = 1;
            case 4
                Y(4,1) = 1;
                Y(8,2) = 1;
        end
        W2(:,:,2*(k1-1)+1) = 1/sqrt(2)*[Y,Y;Y,-Y];
        W2(:,:,2*(k1-1)+2) = 1/sqrt(2)*[Y,Y;1i*Y,-1i*Y];
    end
    W1_back{4} = W1;
    W2_back{4} = W2;
    %% Rank 5
    X = cell(4,1);
    X{1} = 1/2*[1,1,1,1;1,1i,-1,-1i;1,-1,1,-1;1,-1i,-1,1i];
    X{2} = diag([1,exp(1i*pi/4),1i,exp(1i*3*pi/4)])*X{1};
    X{3} = diag([1,exp(1i*pi/8),exp(1i*2*pi/8),exp(1i*3*pi/8)])*X{1};
    X{4} = diag([1,exp(1i*3*pi/8),exp(1i*6*pi/8),exp(1i*9*pi/8)])*X{1};
    W1 = zeros(8,8,4);
    for k = 1:4
        W1(:,:,k) = [X{k},zeros(4,4);zeros(4,4),X{k}];
    end
    e1 = zeros(4,1);
    e1(1) = 1;
    e2 = zeros(4,1);
    e2(2) = 1;
    e3 = zeros(4,1);
    e3(3) = 1;
    e4 = zeros(4,1);
    e4(4) = 1;
    W2 = 1/sqrt(2)*[e1,e1,e2,e2,e3;e1,-e1,e2,-e2,e3];
    W1_back{5} = W1;
    W2_back{5} = W2;
    %% Rank 6
    W1_back{6} = W1;
    W2 = 1/sqrt(2)*[e1,e1,e2,e2,e3,e3;e1,-e1,e2,-e2,e3,-e3];
    W2_back{6} = W2;
    %% Rank 7
    W1_back{7} = W1;
    W2 = 1/sqrt(2)*[e1,e1,e2,e2,e3,e3,e4;e1,-e1,e2,-e2,e3,-e3,e4];
    W2_back{7} = W2;
    %% Rank 8
    W1_back{8} = W1(:,:,1);
    W2 = 1/sqrt(2)*[e1,e1,e2,e2,e3,e3,e4,e4;e1,-e1,e2,-e2,e3,-e3,e4,-e4];
    W2_back{8} = W2;
    
    index = 1;
    for iLayer=1:size(W1_back,1)
        W1 = W1_back{iLayer};
        W2 = W2_back{iLayer};
        for i=1:size(W1,3)
            for j=1:size(W2,3)
                W{index,1} = W1(:,:,i)*W2(:,:,j);
                index = index+1;
            end % end cophasing
        end % end beam selection    
    end % end Layers
end

