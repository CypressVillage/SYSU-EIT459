classdef GFDM
    properties (SetAccess = private)
        SymbolMapping
        Method
        Nr
        PrototypeFilter
        Receiver
    end
    methods
        function obj = GFDM(varargin)
            %% Initialize parameters, set default values
            if numel(varargin)>=5
                obj.Nr.Subcarriers = varargin{1};
                obj.Nr.ActiveSubcarriers = varargin{2};
                obj.Nr.GFDMSymbols = varargin{3};
                obj.Nr.BlockSize = varargin{4};
                obj.Nr.NumberOfBlocks = varargin{5};
            else
                obj.Nr.Subcarriers = 72;
                obj.Nr.ActiveSubcarriers = obj.Nr.Subcarriers;
                obj.Nr.GFDMSymbols = 15;
                obj.Nr.BlockSize = obj.Nr.Subcarriers*obj.Nr.GFDMSymbols;
                obj.Nr.NumberOfBlocks = 1;
            end
            
            if numel(varargin)>=6
                obj.Method = varargin{6};
            else
                obj.Method = 'TTI';
            end
            if numel(varargin)>=7
                obj.PrototypeFilter.PrototypeFilter = varargin{7};
            else
                obj.PrototypeFilter.PrototypeFilter = 'RC';
            end
            
            if numel(varargin)>=8
                obj.PrototypeFilter.RollOff = varargin{8};
            else
                obj.PrototypeFilter.RollOff = 0.1;
            end
            
            if numel(varargin)>=9
                obj.Receiver = varargin{9};
            else
                obj.Receiver = 'ZeroForcing';
            end
            
            switch obj.Method
                case 'TTI'
                    obj.Nr.Subcarriers = 72;
                    obj.Nr.GFDMSymbols    = 15;
                    obj.Nr.BlockSize =  obj.Nr.Subcarriers* obj.Nr.GFDMSymbols;
            end
            
            switch obj.PrototypeFilter.PrototypeFilter
                case 'RRC'
                    g = rrc(obj.Nr.Subcarriers, obj.Nr.GFDMSymbols, obj.PrototypeFilter.RollOff);
                case 'RC'
                    g = rc(obj.Nr.Subcarriers, obj.Nr.GFDMSymbols, obj.PrototypeFilter.RollOff);
                case 'Dirichlet'
                    g = dirichlet(obj.Nr.Subcarriers, obj.Nr.GFDMSymbols);
            end
            
            
        end
         function A_pseudo = Pseudo_matrix(varargin)
             obj = varargin{1};
             A = MatrixA(obj);
             switch obj.Receiver
                 case 'Matched'
                 case 'ZeroForcing'
                 A_pseudo = (A'*A)^(-1)*A';
             end
             
         end
        function PulseShapingMatrix = MatrixG
            g = rrc(obj.Nr.Subcarriers, obj.Nr.GFDMSymbols, obj.PrototypeFilter.RollOff);
            tran_ir = g';
            c = [tran_ir(1) fliplr(tran_ir(end-length(tran_ir)+2:end))];
            G_tx = toeplitz(tran_ir,c);
        end
        
        % taken over from the code of TUDresden
        function A = MatrixA(varargin)
            obj = varargin{1};
            g = rrc(obj.Nr.Subcarriers, obj.Nr.GFDMSymbols, obj.PrototypeFilter.RollOff);
            A = zeros(obj.Nr.Subcarriers*obj.Nr.GFDMSymbols, obj.Nr.Subcarriers*obj.Nr.GFDMSymbols);
            n = 0:obj.Nr.GFDMSymbols*obj.Nr.Subcarriers-1;
            n=n';
            w = exp(1j*2*pi/obj.Nr.Subcarriers);
            
            for k=0:obj.Nr.Subcarriers-1
                for m=0:obj.Nr.GFDMSymbols-1
                    A(:,m*obj.Nr.Subcarriers+k+1) = circshift(g, m*obj.Nr.Subcarriers) .* w.^(k*n);
                end
            end
        end
    end
end

% this part is taken over from TU Dresden
function g = rrc(K, M, a)
% RRC - Return Root Raised Cosine filter (time domain)
t = linspace(-M/2, M/2, M*K+1);
t = t(1:end-1); t = t';
g = (sin(pi*t*(1-a))+4*a.*t.*cos(pi*t*(1+a)))./(pi.*t.*(1-(4*a*t).^2));
g(find(t==0)) = 1-a+4*a/pi;
g(find(abs(t) == 1/(4*a))) = a/sqrt(2)*((1+2/pi)*sin(pi/(4*a))+(1-2/pi)*cos(pi/(4*a)));

g = fftshift(g);
g = g / sqrt(sum(g.*g));
end

function g = rc(K, M, a)
% RC - Return Raised Cosine filter (time domain)
t = linspace(-M/2, M/2, M*K+1);
t = t(1:end-1); t = t';
g = (sinc(t) .* cos(pi*a*t) ./ (1-4*a*a*t.*t));
g = fftshift(g);
g(K+1:K:end) = 0;

g = g / sqrt(sum(g.*g));
end

function g = dirichlet(M, K)
% DIRICHLET - Return Dirichlet Kernel

o = ones(1, M);
z = zeros(1, K*M-M);

G = [o z];
G = circshift(G, [0, -floor(M/2)]);

g = ifft(G);
g = g / sqrt(sum(abs(g).^2));
g = g';
end




