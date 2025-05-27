function subSizes = splitIntoSubarrays(fullSize,numPorts)
% splitIntoSubarrays   Partition a full antenna-array Size vector into numPorts identical sub-sizes
%   fullSize = [nRows, nCols, nPol, nRowPanels, nColPanels]
%   numPorts  = number of CSI-RS ports (i.e. sub-arrays)
%
%   subSizes is an numPorts×5 matrix, each row a [nRows nCols nPol subRowPanels subColPanels].

    % unpack
    n   = fullSize(1);
    m   = fullSize(2);
    p   = fullSize(3);
    ns  = fullSize(4);
    ms  = fullSize(5);

    % how many panels per subarray?
    panelsPerSub = (ns*ms)/numPorts;
    assert(mod(panelsPerSub,1)==0, ...
        'Cannot split %d panels into %d subarrays', ns*ms, numPorts);

    % by default, split along the row-panel direction if possible
    if mod(ns, panelsPerSub)==0
      nsSub = panelsPerSub;
      msSub = 1;
    elseif mod(ms, panelsPerSub)==0
      nsSub = 1;
      msSub = panelsPerSub;
    else
      error('Unable to partition %d×%d panel grid into %d chunks', ns, ms, numPorts)
    end

    subSizes = repmat([n m p nsSub msSub], numPorts, 1);
end