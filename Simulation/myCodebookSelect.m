
% function W = myCodebookSelect(pmi, nTx, nLayers)
% % myCodebookSelect - Custom simplified Type 1 Single Panel codebook precoder
% % Supports 4 or 8 Tx antennas and 1, 2, 4, or 8 layers
% %
% % Parameters:
% %   pmi - Precoding Matrix Index (PMI)
% %   nTx - Number of Tx antennas (4 or 8)
% %   nLayers - Number of transmission layers (1, 2, 4, or 8)
% 
% % Validate inputs
% if ~(nTx == 4 || nTx == 8)
%     error('Only 4 or 8 Tx antennas are supported.');
% end
% 
% % Identity-style mappings (you can expand these with actual 3GPP codebook options)
% if nLayers == 1
%     W = zeros(1, nTx);
%     W(1) = 1;
% 
% elseif nLayers == 2
%     if nTx == 4
%         codebook = {
%             [1 0 0 0; 0 1 0 0], ...
%             [0 1 0 0; 0 0 1 0], ...
%             [0 0 1 0; 0 0 0 1]
%         };
%     else % 8 Tx
%         codebook = {
%             [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0], ...
%             [1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0], ...
%             [1 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0]
%         };
%     end
%     idx = mod(pmi, length(codebook)) + 1;
%     W = codebook{idx};
% 
% elseif nLayers == 4
%     if nTx == 4
%         W = eye(4);
%     else % 8 Tx
%         codebook = {
%             [eye(4), zeros(4,4)], ...
%             [zeros(4,4), eye(4)]
%         };
%         idx = mod(pmi, length(codebook)) + 1;
%         W = codebook{idx};
%     end
% 
% elseif nLayers == 8
%     if nTx == 8
%         W = eye(8);
%     else
%         error('8 layers require 8 Tx antennas');
%     end
% 
% else
%     error('Unsupported layer configuration. Only 1, 2, 4, or 8 layers supported.');
% end
% end

function W = myCodebookSelect(pmi, nTx, nLayers)
% myCodebookSelect - 3GPP TS 38.214 Type-1 Single-Panel PMI codebook precoder
%   W = myCodebookSelect(pmi, nTx, nLayers)
%     returns the nTx×nLayers precoder for the given zero-based PMI index,
%     number of CSI-RS ports nTx∈{4,8}, and number of layers nLayers∈{1,2,4,8}.

  %% validate inputs
  if ~isscalar(pmi) || pmi<0 || floor(pmi)~=pmi
    error('PMI must be a non‐negative integer scalar.');
  end
  if ~ismember(nTx,[4 8])
    error('nTx must be 4 or 8.');
  end
  if ~ismember(nLayers,[1 2 4 8])
    error('nLayers must be 1, 2, 4, or 8.');
  end
  if nLayers>nTx
    error('Number of layers cannot exceed number of ports.');
  end

  %% pick the (N1,N2,O1,O2) per TS 38.214 Tbl 5.2.2.2.1-2
  switch nTx
    case 4,  N1=2; N2=1; O1=4; O2=1;
    case 8,  N1=4; N2=1; O1=4; O2=1;
  end

  %% generate the full codebook for this (nTx,nLayers)
  switch nLayers
    case 1, Wset = gen1layer(nTx,N1,N2,O1,O2);
    case 2, Wset = gen2layer(nTx,N1,N2,O1,O2);
    case 4, Wset = gen4layer(nTx,N1,N2,O1,O2);
    case 8, Wset = gen8layer(nTx,N1,N2,O1,O2);
  end

  %% extract the requested PMI entry
  idx = pmi + 1;  % MATLAB 1-based
  if idx > size(Wset,3)
    error('PMI index %d out of range for %d-layer, %d-port codebook.',...
          pmi,nLayers,nTx);
  end
  W = Wset(:,:,idx);
end


%% ───  Local helper functions ──────────────────────────────────────────────

function W = gen1layer(P,N1,N2,O1,O2)
  % Table 5.2.2.2.1-5
  M1 = N1*O1;  M2 = N2*O2;
  cnt = 0;
  for i11 = 0:M1-1
    for i12 = 0:M2-1
      for i2 = 0:2
        cnt = cnt+1;
        w = zeros(P,1);
        for m = 1:P
          % horizontal + vertical steering
          alpha = (m-1)*i11/M1;
          beta  = ceil(m/N1-1)*i12/M2;
          w(m) = exp(-1j*2*pi*(alpha + beta));
        end
        W(:,:,cnt) = w / sqrt(P);
      end
    end
  end
end

function W = gen2layer(P,N1,N2,O1,O2)
  % Table 5.2.2.2.1-6
  M1 = N1*O1;  M2 = N2*O2;
  cnt = 0;
  for i11 = 0:M1-1
    for i12 = 0:M2-1
      for i2 = 0:1
        cnt = cnt+1;
        w = zeros(P,2);
        for m = 1:P
          alpha = (m-1)*i11/M1 + ceil(m/N1-1)*i12/M2;
          w(m,1) = exp(-1j*2*pi*alpha);
          w(m,2) = exp(-1j*2*pi*(alpha + i2/M1));
        end
        W(:,:,cnt) = w / sqrt(2*P);
      end
    end
  end
end

function W = gen4layer(P,N1,N2,O1,O2)
  % Table 5.2.2.2.1-8
  M1 = N1*O1;  M2 = N2*O2;
  cnt = 0;
  for i11 = 0:M1-1
    for i12 = 0:M2-1
      for i13 = 0:1
        for i2 = 0:1
          cnt = cnt+1;
          w = zeros(P,4);
          for m = 1:P
            alpha = (m-1)*i11/M1 + ceil(m/N1-1)*i12/M2;
            w(m,1) = exp(-1j*2*pi* alpha           );
            w(m,2) = exp(-1j*2*pi*(alpha + i2/M1)  );
            w(m,3) = exp(-1j*2*pi*(alpha + i13/M1) );
            w(m,4) = exp(-1j*2*pi*(alpha + (i13+i2)/M1));
          end
          W(:,:,cnt) = w / sqrt(4*P);
        end
      end
    end
  end
end

function W = gen8layer(P,N1,N2,O1,O2)
  % Table 5.2.2.2.1-12  (only valid for P=8,N1=4,N2=1,O1=4,O2=1)
  M1 = N1*O1;  M2 = N2*O2;
  cnt = 0;
  for i11 = 0:(M1/2)-1
    for i12 = 0:M2-1
      for i13 = 0:3
        for i2 = 0:1
          cnt = cnt+1;
          w = zeros(P,8);
          for m = 1:P
            alpha = (m-1)*i11/(M1/2) + ceil(m/N1-1)*i12/M2;
            % first 4 layers
            for k = 0:3
              w(m,k+1)   = exp(-1j*2*pi*(alpha + k*i2/(M1/2)));
              w(m,k+5) = exp(-1j*2*pi*(alpha + i13/(M1/2) + k*i2/(M1/2)));
            end
          end
          W(:,:,cnt) = w / sqrt(8*P);
        end
      end
    end
  end
end