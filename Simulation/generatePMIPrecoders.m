function Wband = generatePMIPrecoders(pmiMode, fixedPMI, nTx, nLayers, nBands)
%GENERATEPMIPRECODERS  Create PMI‐based precoders for wideband or subband.
%   Wband = generatePMIPrecoders(pmiMode, fixedPMI, nTx, nLayers, nBands)
%   returns a [nTx x nLayers × nBands] array of precoding matrices.
%
%   pmiMode   — 'wideband' or 'subband'
%   fixedPMI  — [] to randomize, or a scalar (wideband) / vector (subband)
%   nTx       — number of CSI-RS ports (4 or 8)
%   nLayers   — number of layers (1,2,4,8)
%   nBands    — number of subbands (for 'subband' mode)

  % Pre-allocate 
  Wband = zeros(nTx, nLayers, nBands);

  switch lower(pmiMode)
    case 'wideband'
      % single PMI for all bands
      if isempty(fixedPMI)
        pmi = randi([0, nPMI(nTx,nLayers)-1]);
      else
        pmi = fixedPMI;
      end
      W = myCodebookSelect(pmi, nTx, nLayers);  % nTx×nLayers
      for b = 1:nBands
        Wband(:,:,b) = W;
      end

    case 'subband'
      % independent PMI per band
      if isempty(fixedPMI)
        for b = 1:nBands
          pmi = randi([0, nPMI(nTx,nLayers)-1]);
          Wband(:,:,b) = myCodebookSelect(pmi, nTx, nLayers);
        end
      else
        assert(numel(fixedPMI)==nBands, ...
          'fixedPMI must be a vector of length nBands for subband mode.');
        for b = 1:nBands
          Wband(:,:,b) = myCodebookSelect(fixedPMI(b), nTx, nLayers);
        end
      end

    otherwise
      error('Unsupported pmiMode: %s', pmiMode);
  end
end