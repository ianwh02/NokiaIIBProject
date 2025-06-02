simParameters = struct();       % Clear simParameters variable to contain all key simulation parameters 
simParameters.NFrames = 1;      % Number of 10 ms frames
simParameters.SNRIn = 25; % SNR range (dB)
simParameters.PerfectChannelEstimator = true;
simParameters.DisplaySimulationInformation = true;
simParameters.DisplayDiagnostics = false;

% SCS carrier parameters
simParameters.Carrier = nrCarrierConfig;         % Carrier resource grid configuration
simParameters.Carrier.NSizeGrid = 52;            % Bandwidth in number of resource blocks
simParameters.Carrier.SubcarrierSpacing = 30;    % 15, 30, 60, 120 (kHz)
simParameters.Carrier.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
simParameters.Carrier.NCellID = 1;               % Cell identity

simParameters.PDSCH = nrPDSCHConfig;      % This PDSCH definition is the basis for all PDSCH transmissions in the simulation

% Define PDSCH time allocation in a slot
simParameters.PDSCH.MappingType = "A"; % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
symAlloc = [2 simParameters.Carrier.SymbolsPerSlot-2];  % Starting symbol and number of symbols of each PDSCH allocation
simParameters.PDSCH.SymbolAllocation = symAlloc;

% Define PDSCH frequency resource allocation per slot to be full grid 
simParameters.PDSCH.PRBSet = 0:simParameters.Carrier.NSizeGrid-1;

% Scrambling identifiers
simParameters.PDSCH.NID = simParameters.Carrier.NCellID;
simParameters.PDSCH.RNTI = 1;

simParameters.PDSCH.NumLayers = 4;

simParameters.PDSCH.DMRS.DMRSTypeAPosition       = 2; % Mapping type A only. First DM-RS symbol position (2,3)
simParameters.PDSCH.DMRS.DMRSLength              = 1; % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
simParameters.PDSCH.DMRS.DMRSAdditionalPosition  = 1; % Additional DM-RS symbol positions (max range 0...3)
simParameters.PDSCH.DMRS.DMRSConfigurationType   = 1; % DM-RS configuration type (1,2)
simParameters.PDSCH.DMRS.NumCDMGroupsWithoutData = 3; % CDM groups without data (1,2,3)

simParameters.PDSCHExtension = struct();
simParameters.PDSCHExtension.PRGBundleSize = 2; % 2, 4, or [] to signify "wideband"
simParameters.PDSCHExtension.MCSTable      = 'Table1'; % 'Table1',...,'Table4'
simParameters.PDSCHExtension.XOverhead     = [ ]; % 0, 6, 12, 18, or [] for automatic selection.
simParameters.PDSCHExtension.NHARQProcesses   = 8;       % e.g. 8 parallel HARQ processes
simParameters.PDSCHExtension.EnableHARQ       = true;    % enable retransmissions
simParameters.PDSCHExtension.EnableCBGTransmission = true; % Enable CBG-based transmission, otherwise TB-based transmission
simParameters.PDSCHExtension.MaxNumCBG = 4;          % Maximum number of CBGs per transport block for each HARQ process in CBG-based transmission
simParameters.PDSCHExtension.RVSequence       = [0 2 3 1]; % standard 4-RV cycle
simParameters.PDSCHExtension.TargetCodeRate = 567/1024;

simParameters.PDSCHExtension.LDPCDecodingAlgorithm = "Layered belief propagation";
simParameters.PDSCHExtension.MaximumLDPCIterationCount = 6;

simParameters.TransmitAntennaArray = struct('Size',[8 4 2 8 1], ... % [8 2 2 8 1] for 4Tx, [8 4 2 8 1] for 8Tx
        'ElementSpacing',[0.5 0.5 4 1], ...
        'PolarizationAngles', [-45 45], ...
        'Orientation', [0; 10; 0], ...
        'Element', '38.901', ...
        'PolarizationModel', 'Model-2'); 
simParameters.ReceiveAntennaArray = struct('Size',[2 1 2 1 1], ... % [2 1 2 1 1] for 4Rx, [2 2 2 1 1] for 8Rx
        'ElementSpacing',[0.5 0.5 1 1], ...
        'PolarizationAngles', [-45 45], ...
        'Orientation', [180; 0; 0], ...
        'Element', 'isotropic', ...
        'PolarizationModel', 'Model-2');

simParameters.TransmitAntennaArray.NumPanels        = 4; % Number of panels in horizontal dimension (Ng)
simParameters.TransmitAntennaArray.PanelDimensions  = [4 1]; % Number of columns and rows in a panel (N1, N2)
simParameters.TransmitAntennaArray.NumPolarizations = 2; % Number of transmit polarizations
simParameters.ReceiveAntennaArray.NumPanels         = 1; % Number of panels in horizontal dimension (Ng)
simParameters.ReceiveAntennaArray.PanelDimensions   = [2 2];                % Number of columns and rows in a panel (N1, N2)
simParameters.ReceiveAntennaArray.NumPolarizations  = 2; % Number of receive polarizations

simParameters.NTxAnts = numAntennaElements(simParameters.TransmitAntennaArray);
simParameters.NRxAnts = numAntennaElements(simParameters.ReceiveAntennaArray);

simParameters.CSIRS = nrCSIRSConfig;
simParameters.CSIRS.CSIRSType = 'nzp'; % 'nzp','zp'
simParameters.CSIRS.RowNumber = 16; % 4 for 4 CSI, 6 for 8 CSI
simParameters.CSIRS.NumRB = simParameters.Carrier.NSizeGrid - simParameters.CSIRS.RBOffset;
simParameters.CSIRS.CSIRSPeriod = 'on';
simParameters.CSIRS.SymbolLocations = [2, 4];
simParameters.CSIRS.SubcarrierLocations = [0,3,6,9]; % 0 or [0,3,6,9]
simParameters.CSIRS.Density = 'one';

disp(['Number of CSI-RS ports: ' num2str(simParameters.CSIRS.NumCSIRSPorts) '.'])
csirsCDMLengths = getCSIRSCDMLengths(simParameters.CSIRS);

% Check that the number of CSI-RS ports and transmit antenna elements match
% and the consistency of multiple CSI-RS resources
validateCSIRSConfig(simParameters.Carrier,simParameters.CSIRS,simParameters.NTxAnts);

simParameters.CSIReportMode = 'RI-PMI-CQI'; % 'RI-PMI-CQI','AI CSI compression','Perfect CSI'

simParameters.CSIReportConfig = struct();
simParameters.CSIReportConfig.Period = [5 0];  % Peridocity and offset of the CSI report in slots

if simParameters.CSIReportMode == "RI-PMI-CQI"  
    
    riRestrict = zeros(1,4);
    riRestrict(simParameters.PDSCH.NumLayers) = 1;

    simParameters.CSIReportConfig.CQITable          = "Table3"; % 'Table1','Table2','Table3'
    simParameters.CSIReportConfig.CQIMode           = 'Wideband'; % 'Wideband','Subband'
    simParameters.CSIReportConfig.PMIMode           = 'Subband'; % 'Wideband','Subband'
    simParameters.CSIReportConfig.CodebookType      = 'Type1MultiPanel'; % 'Type1SinglePanel','Type1MultiPanel','Type2'
    simParameters.CSIReportConfig.SubbandSize       = 4; % Subband size in RB (4,8,16,32)
    simParameters.CSIReportConfig.CodebookMode      = 1; % 1,2
    simParameters.CSIReportConfig.RIRestriction     = riRestrict;                   % Empty for no rank restriction
    simParameters.CSIReportConfig.NumberOfBeams     = 2; % 2,3,4. Only for Type II codebooks
    simParameters.CSIReportConfig.PhaseAlphabetSize = 8; % 4,8. Only for Type II codebooks
    simParameters.CSIReportConfig.SubbandAmplitude  = true;                  % true/false. Only for Type II codebooks
    simParameters.CSIReportConfig.NStartBWP         = [];                   % Empty to signal the entire carrier
    simParameters.CSIReportConfig.NSizeBWP          = [];                   % Empty to signal the entire carrier
    
    simParameters.CSIReportConfig.PMIModeOverride = 'random';  % 'best','random','fixed'
    simParameters.CSIReportConfig.FixedPMI       = 3;       % zero‐based PMI if you choose 'fixed'

    % Configure the CSI report with the antenna panel dimensions specified
    simParameters.CSIReportConfig.PanelDimensions = getCSIReportPanelDimensions(simParameters.TransmitAntennaArray,simParameters.CSIReportConfig.CodebookType);
    
    % Adjust the RIRestriction according to the provided DM-RS
    % configuration. If the number of DM-RS antenna ports in a given
    % configuration is less than the possible number of ranks for a CSI
    % report, restrict those ranks through CSIReportConfig
    simParameters.CSIReportConfig.RIRestriction = updateRankRestriction(simParameters.PDSCH.DMRS,simParameters.CSIReportConfig);
else % AI CSI compression
    
    % Specify the file name of the AI neural network
    simParameters.AINetworkFilename = 'csiTrainedNetwork.mat'; 

end

simParameters.UEProcessingDelay = 7;
simParameters.BSProcessingDelay = 1;

simParameters.DelayProfile = 'CDL-C';   % 'CDL-A',...,'CDL-E','TDL-A',...,'TDL-E'
simParameters.DelaySpread = 300e-9;     % s
simParameters.MaximumDopplerShift = 10;  % Hz

[channel, simParameters] = createChannel(simParameters);
simParameters.Channel = channel;

disp("CDL Tx array size: " + mat2str(simParameters.Channel.TransmitAntennaArray.Size));
nTxAnts = prod(simParameters.Channel.TransmitAntennaArray.Size);
disp("Number of physical Tx elements: " + nTxAnts);
disp("WtxVirtual size: " + mat2str(size(simParameters.WtxVirtual)));

% Create DL-SCH coder / decoder with HARQ
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate         = simParameters.PDSCHExtension.TargetCodeRate;

decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate         = simParameters.PDSCHExtension.TargetCodeRate;
decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

errorBlocksPar = zeros(numel(simParameters.SNRIn),1);
totalBlocksPar = zeros(numel(simParameters.SNRIn),1);
simThroughputPar = zeros(numel(simParameters.SNRIn),1);
maxThroughputPar = zeros(numel(simParameters.SNRIn),1);
CSIReportPar = cell(numel(simParameters.SNRIn),1);
blerPar = zeros(numel(simParameters.SNRIn),1);

for snrIdx = 1:numel(simParameters.SNRIn)
%parfor snrIdx = 1:numel(simParameters.SNRIn)
% To reduce the total simulation time, you can execute this loop in
% parallel by using the Parallel Computing Toolbox. Comment out the 'for'
% statement and uncomment the 'parfor' statement.
    
    % Reset the random number generator for repeatability
    rng(snrIdx + 1000, "twister")

    % Display simulation information at this SNR point
    displaySNRPointProgress(simParameters,snrIdx);

    % Take full copies of the simulation-level parameter structures so that
    % they are not PCT broadcast variables when using parfor
    simParamLocal = simParameters;

    % Set up the transmitter, propagation channel, receiver, and CSI
    % feedback configuration parameters.
    [carrier,encodeDLSCH,pdsch,pdschextra,csirs,wtx] = setupTransmitter(simParamLocal);
    [channel,maxChDelay,N0] = setupChannel(simParamLocal,snrIdx);
    [decodeDLSCH,pathFilters,timingOffset] = setupReceiver(simParamLocal);
    csiFeedbackOpts = getCSIFeedbackOptions(simParamLocal,snrIdx);

    % Obtain an initial CSI report based on perfect channel estimates that
    % the Tx can use to adapt the transmission parameters.
    [csiReports, Hest] = initialCSIReport(simParamLocal,snrIdx,carrier,csirs,channel,csiFeedbackOpts);
    csiAvailableSlots = 0;
    
    if simParameters.CSIReportMode == "RI-PMI-CQI" 

        mode      = simParameters.CSIReportConfig.PMIModeOverride;
        reportCfg = simParameters.CSIReportConfig;
        nLayers   = pdsch.NumLayers;
    
        switch lower(mode)
          case 'best'
            % do nothing, keep the UE‐chosen W in csiReports.W
        
          case 'random'
            % replace with truly random PMI
            [~,info]      = hDLPMIRandom(carrier, csirs, reportCfg, nLayers, Hest);
            csiReports.W   = info.W;
    
          otherwise
            error('Unknown PMI override mode "%s".',mode);
        end
    end

    % Total number of slots in the simulation period
    NSlots = simParamLocal.NFrames * carrier.SlotsPerFrame;
    
    % Set up HARQ sequencing
    harqSeq   = 0:(simParameters.PDSCHExtension.NHARQProcesses-1);
    rvSeq     = simParameters.PDSCHExtension.RVSequence;
    harqEntity = HARQEntity(harqSeq, rvSeq, simParameters.PDSCH.NumCodewords);

    % Loop over the entire waveform length
    for nslot = 0:NSlots-1

        % Update the carrier slot numbers for new slot
        carrier.NSlot = nslot;

        % Use new CSI report to configure the number of layers and
        % modulation of the PDSCH and target code rate of the DL-SCH if
        % there is a new report available.
        [isNewReport,repIdx] = ismember(nslot,csiAvailableSlots);
        if isNewReport
            [pdsch.Modulation,pdschextra.TargetCodeRate,wtx] = hCSIDecode(carrier,pdsch,pdschextra,csiReports(repIdx),csiFeedbackOpts);
            pdsch.NumLayers = size(wtx,1);
            pdsch.Modulation = '64QAM';
            if pdsch.NumCodewords > 1
                pdschextra.TargetCodeRate = [567/1024, 567/1024];
            else 
                pdschextra.TargetCodeRate = 567/1024;
            end
            encodeDLSCH.TargetCodeRate = pdschextra.TargetCodeRate;
        end
        % Create an OFDM resource grid for a slot
        dlGrid = nrResourceGrid(carrier,csirs.NumCSIRSPorts);

        % CSI-RS mapping to the slot resource grid
        [csirsInd,csirsInfo] = nrCSIRSIndices(carrier,csirs);
        csirsSym = nrCSIRS(carrier,csirs);

        dlGrid(csirsInd) = csirsSym;
        csirsTransmission = ~isempty(csirsInd);

        % PDSCH reserved REs for CSI-RS
        pdsch.ReservedRE = csirsInd-1; % 0-based indices
        
        % PDSCH generation
        % Calculate the transport block sizes for the transmission in the slot
        [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschextra.TargetCodeRate,pdschextra.XOverhead);

        % ─── HARQ-controlled TB generation & encoding ───
        for cwIdx = 1:pdsch.NumCodewords
            procID = harqEntity.HARQProcessID;       % current HARQ process
            rv     = harqEntity.RedundancyVersion(cwIdx);
        
            % New data only if this is a fresh transmission
            if harqEntity.NewData(cwIdx)
                trBlk = randi([0 1], trBlkSizes(cwIdx), 1);
                setTransportBlock(encodeDLSCH, trBlk, cwIdx-1, procID);
            end
        
            % If this process has timed out, clear its soft buffer
            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCH, cwIdx-1, procID);
            end
        end
        
        % Now encode using the per-CW RVs and process ID
        RVs = harqEntity.RedundancyVersion;  % 1×NumCodewords
        codedTrBlocks = encodeDLSCH(pdsch.Modulation, pdsch.NumLayers, pdschIndicesInfo.G, RVs);

        % PDSCH modulation and precoding
        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlocks);
        [pdschAntSymbols,pdschAntIndices] = nrPDSCHPrecode(carrier,pdschSymbols,pdschIndices,wtx);
        dlGrid(pdschAntIndices) = pdschAntSymbols;        

        % PDSCH DM-RS precoding and mapping
        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
        [dmrsAntSymbols,dmrsAntIndices] = nrPDSCHPrecode(carrier,dmrsSymbols,dmrsIndices,wtx);
        dlGrid(dmrsAntIndices) = dmrsAntSymbols;

        % OFDM modulation
        txWaveform = nrOFDMModulate(carrier, dlGrid);  % [N × 8], where 8 = virtual ports
        
        txWaveform = txWaveform * simParameters.WtxVirtual.';

        % Pad to flush channel memory
        txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform,2))];
        
        % Pass through 512-element physical channel
        [rxWaveform, pathGains, sampleTimes] = channel(txWaveform);
        
        % Add AWGN to the received time-domain waveform
        noise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
        rxWaveform = rxWaveform + noise;

        if (simParamLocal.PerfectChannelEstimator)
            % Perfect synchronization. Use information provided by the
            % channel to find the strongest multipath component
            pathFilters = getPathFilters(channel);
            [timingOffset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
        else
            % Practical synchronization. Correlate the received waveform
            % with the PDSCH DM-RS to obtain the timing offset and
            % correlation magnitude. The receiver updates the timing offset
            % only when the correlation magnitude is high.
            [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols); 
            timingOffset = hSkipWeakTimingOffset(timingOffset,t,mag);
            % Display a warning if the estimated timing offset exceeds the
            % maximum channel delay
            if timingOffset > maxChDelay
                warning(['Estimated timing offset (%d) is greater than the maximum channel delay (%d).' ...
                    ' This will result in a decoding failure. This may be caused by low SNR,' ...
                    ' or not enough DM-RS symbols to synchronize successfully.'],timingOffset,maxChDelay);
            end
        end
        rxWaveform = rxWaveform(1+timingOffset:end,:);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
        [K,L,R] = size(rxGrid);
        if (L < carrier.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,carrier.SymbolsPerSlot-L,R));
        end

        if (simParamLocal.PerfectChannelEstimator)
            % Perfect channel estimation, using the value of the path gains
            % provided by the channel. This channel estimate does not
            % include the effect of transmitter precoding
            Hest = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,timingOffset,sampleTimes);
            % Apply subarray virtualization: reduce 512 → 8 ports
            Hest = virtualizeChannel(Hest, simParameters.WtxVirtual);

            % Get perfect noise estimate (from the noise realization)
            noiseGrid = nrOFDMDemodulate(carrier,noise(1+timingOffset:end ,:));
            noiseEst = var(noiseGrid(:));

            % Get PDSCH resource elements from the received grid and 
            % channel estimate
            [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,Hest);

            % Apply precoding to channel estimate
            pdschHest = nrPDSCHPrecode(carrier,pdschHest,pdschHestIndices,permute(wtx,[2 1 3]));
        else
            % Practical channel estimation between the received grid and
            % each transmission layer, using the PDSCH DM-RS for each
            % layer. This channel estimate includes the effect of
            % transmitter precoding
            [Hest,noiseEst] = hSubbandChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,pdschextra.PRGBundleSize,'CDMLengths',pdsch.DMRS.CDMLengths);

            % Average noise estimate across PRGs and layers
            noiseEst = mean(noiseEst,'all');
            
            % Get PDSCH resource elements from the received grid and
            % channel estimate
            [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,Hest);
        end

        % Equalization
        [pdschEq,eqCSIScaling] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        % Decode PDSCH physical channel
        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);
        
        % Display EVM per layer, per slot and per RB
        if (simParamLocal.DisplayDiagnostics)
            plotLayerEVM(NSlots,nslot,pdsch,size(dlGrid),pdschIndices,pdschSymbols,pdschEq);
        end
        
        % Scale LLRs
        eqCSIScaling = nrLayerDemap(eqCSIScaling); % CSI scaling layer demapping
        for cwIdx = 1:pdsch.NumCodewords
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx});        % bits per symbol
            eqCSIScaling{cwIdx} = repmat(eqCSIScaling{cwIdx}.',Qm,1);      % expand by each bit per symbol
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* eqCSIScaling{cwIdx}(:); % scale LLRs
        end
        
        % Decode DL-SCH
        decodeDLSCH.TransportBlockLength = trBlkSizes;
        [rxBits, blkerr] = decodeDLSCH(dlschLLRs, pdsch.Modulation, pdsch.NumLayers, harqEntity.RedundancyVersion);
        
        % Update HARQ state
        status = updateAndAdvance(harqEntity, blkerr, trBlkSizes, pdschIndicesInfo.G);
        if simParamLocal.DisplaySimulationInformation
            fprintf(' [%s]', status);
        end

        % Treat any codeword error as a block error:
        errorBlocksPar(snrIdx) = errorBlocksPar(snrIdx) + any(blkerr);
        totalBlocksPar(snrIdx) = totalBlocksPar(snrIdx) + 1;

        % Store values to calculate throughput
        simThroughputPar(snrIdx) = simThroughputPar(snrIdx) + sum(~blkerr .* trBlkSizes);
        maxThroughputPar(snrIdx) = maxThroughputPar(snrIdx) + sum(trBlkSizes);

        % CSI measurements and encoding 
        if csirsTransmission
            
            if (~simParamLocal.PerfectChannelEstimator)
                % Consider only the NZP-CSI-RS symbols and indices for CSI-RS based
                % channel estimation
                nzpind = (csirsSym ~= 0);
                
                % Calculate practical channel estimate based on CSI-RS. Use
                % a time-averaging window that covers all of the
                % transmitted CSI-RS symbols.
                [Hest,noiseEst] = nrChannelEstimate(carrier,rxGrid, ...
                    csirsInd(nzpind),csirsSym(nzpind),'CDMLengths',csirsCDMLengths);

            end

            % Generate CSI report. Store the report for use at the
            % transmitter. The CSI feedback is subject to a delay that
            % depends on the CSI report periodicity and the UE processing
            % delay. The slot in which the CSI is available to use at the
            % transmitter depends on the BS processing delay as well.
            rxCSIReport = hCSIEncode(carrier,csirs,Hest,noiseEst,csiFeedbackOpts);
            csiFeedbackSlot = nextCSISlot(csiFeedbackOpts.CSIReportPeriod,1+nslot+simParamLocal.UEProcessingDelay);
            csiAvailableSlots(end+1) = 1+csiFeedbackSlot+simParamLocal.BSProcessingDelay; %#ok<SAGROW>
            csiReports(end+1) = rxCSIReport; %#ok<SAGROW>
        end

        % Print slot-wise information
        if (simParamLocal.DisplaySimulationInformation)
            printSlotInfo(NSlots,carrier,pdsch,pdschextra,blkerr,trBlkSizes./pdschIndicesInfo.G,csirsTransmission,csiReports,repIdx)
        end
    end

    % Store CSI report for each SNR point
    CSIReportPar{snrIdx} = csiReports;

    % After finishing all slots for this SNR
    blerPar(snrIdx) = errorBlocksPar(snrIdx) / totalBlocksPar(snrIdx);
    
    % Display the results dynamically in the command window
    if (simParamLocal.DisplaySimulationInformation)
        fprintf('\n');
    end
    fprintf('\nThroughput(Mbps) for %d frame(s) = %.4f\n',simParamLocal.NFrames,1e-6*simThroughputPar(snrIdx)/(simParamLocal.NFrames*10e-3));
    fprintf('Throughput(%%) for %d frame(s) = %.4f\n',simParamLocal.NFrames,simThroughputPar(snrIdx)*100/maxThroughputPar(snrIdx));
    fprintf('BLER for %d frame(s) = %.4f\n', simParameters.NFrames, blerPar(snrIdx));

end

bler = blerPar;
maxThroughput = maxThroughputPar;
simThroughput = simThroughputPar;
CSIReport = CSIReportPar;

figure;
plot(simParameters.SNRIn,simThroughput*100./maxThroughput,'o-.')

xlabel('SNR (dB)'); ylabel('Throughput (%)'); grid on;
title(sprintf('%s (%dx%d) / NRB=%d / SCS=%dkHz / CSI: %s', ...
              simParameters.DelayProfile,simParameters.NTxAnts,simParameters.NRxAnts, ...
              simParameters.Carrier.NSizeGrid,simParameters.Carrier.SubcarrierSpacing,...
              char(simParameters.CSIReportMode)));
if simParameters.CSIReportMode == "RI-PMI-CQI"
    perc = 90;
    plotCQI(simParameters,CSIReport,perc)
    simResults.CSIReport = CSIReport;
end

% Bundle key parameters and results into a combined structure for recording
simResults.simParameters = simParameters;
simResults.simThroughput = simThroughput;
simResults.maxThroughput = maxThroughput;
simResults.bler         = bler;

resultsFolder = './results/';
filename = sprintf('simResult_1Y8Txsubarray_%s_%dTx_%dRx_%dLayer_%dHz.mat', ...
    simParameters.DelayProfile, simParameters.NTxAnts, simParameters.NRxAnts, simParameters.PDSCH.NumLayers, simParameters.MaximumDopplerShift);
fullPath = fullfile(resultsFolder, filename);

if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

save(fullPath, 'simResults', '-v7.3');

function [carrier,eDLSCH,pdsch,pdschextra,csirs,wtx] = setupTransmitter(simParameters)
% Extract channel and signal-level parameters, create DL-SCH encoder, and
% initialize MIMO precoding matrix.

    carrier = simParameters.Carrier;
    pdsch = simParameters.PDSCH;
    pdschextra = simParameters.PDSCHExtension;
    csirs = simParameters.CSIRS;

    % Select XOverhead for TBS calculation if required
    if isempty(pdschextra.XOverhead)
        pdschextra.XOverhead = getXOverhead(carrier,csirs);
    end
    
    % Create DL-SCH encoder system object to perform transport channel
    % encoding
    eDLSCH = nrDLSCH;
    
    % Normalize columns
    wtx = 1;
    
end

function [decodeDLSCH,pathFilters,timingOffset] = setupReceiver(simParameters)
% Create and configure DL-SCH decoder and initialize receiver parameters

    % Create DL-SCH decoder system object to perform transport channel
    % decoding
    decodeDLSCH = nrDLSCHDecoder;
    decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
    decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

    % Initialize channel path filters and timing offset. Timing offset is
    % updated in every slot for perfect synchronization and when the
    % correlation is strong for practical synchronization
    pathFilters = [];
    timingOffset = 0;

end

function [channel, simParameters] = createChannel(simParameters)
% Create and configure the propagation channel. If the number of antennas
% is 1, configure only 1 polarization, otherwise configure 2 polarizations.

    % Number of antenna elements and polarizations
    nTxAnts = simParameters.NTxAnts;
    numTxPol = simParameters.TransmitAntennaArray.NumPolarizations;
    nRxAnts = simParameters.NRxAnts;
    numRxPol = simParameters.ReceiveAntennaArray.NumPolarizations;
    
    if contains(simParameters.DelayProfile,'CDL')

        % Create CDL channel
        channel = nrCDLChannel;

        % Tx antenna array configuration in CDL channel. The number of antenna
        % elements depends on the panel dimensions. The size of the antenna
        % array is [M,N,P,Mg,Ng]. M and N are the number of rows and columns in
        % the antenna array. P is the number of polarizations (1 or 2). Mg and
        % Ng are the number of row and column array panels respectively. Note
        % that N1 and N2 in the panel dimensions follow a different convention
        % and denote the number of columns and rows, respectively.

        channel.TransmitAntennaArray.Size = simParameters.TransmitAntennaArray.Size;
        channel.TransmitAntennaArray.ElementSpacing = simParameters.TransmitAntennaArray.ElementSpacing;
        channel.TransmitAntennaArray.PolarizationAngles = simParameters.TransmitAntennaArray.PolarizationAngles;
        channel.TransmitAntennaArray.Orientation = simParameters.TransmitAntennaArray.Orientation;
        channel.TransmitAntennaArray.Element = simParameters.TransmitAntennaArray.Element;
        channel.TransmitAntennaArray.PolarizationModel = simParameters.TransmitAntennaArray.PolarizationModel;

        % Extract array dimensions
        antSize = simParameters.TransmitAntennaArray.Size;  % [M N P Mg Ng]
        numElements = prod(antSize);
        
        % Full 5D grid to get 512 positions
        M = antSize(1); N = antSize(2); P = antSize(3);
        Mg = antSize(4); Ng = antSize(5);
        dx = simParameters.TransmitAntennaArray.ElementSpacing(1);
        dy = simParameters.TransmitAntennaArray.ElementSpacing(2);
        dz = 1;  % Assume 1 unit vertical for polarization or stacking
        
        elementPos = zeros(3, prod(antSize));  % Preallocate [3 × 512]
        idx = 1;
        
        for ng = 0:Ng-1
            for mg = 0:Mg-1
                for p = 0:P-1
                    for n = 0:N-1
                        for m = 0:M-1
                            x = (n + N*ng) * dx;
                            y = (m + M*mg) * dy;
                            z = p * dz;
                            elementPos(:,idx) = [x; y; z];
                            idx = idx + 1;
                        end
                    end
                end
            end
        end
        
        % Manual subarray mapping
        numVirtualPorts = simParameters.CSIRS.NumCSIRSPorts;
        elementIndices = (1:numElements).';
        subarrayMap = mat2cell( ...
            elementIndices, ...
            repmat(floor(numElements/numVirtualPorts), numVirtualPorts, 1), ...
            1);
        
        remaining = numElements - sum(cellfun(@numel, subarrayMap));
        if remaining > 0
            subarrayMap{end} = [subarrayMap{end}; (numElements - remaining + 1 : numElements)'];
        end
        
        % Steering vectors using manual position-based approach
        steeringVecs = zeros(numElements, numVirtualPorts);
        
        azAngles = linspace(-45, 45, 4);  % Horizontal directions
        elAngles = linspace(-5, -25, 2);  % Elevation: shallow to steep downward
        [AZ, EL] = meshgrid(azAngles, elAngles);
        steeringDir = [AZ(:)'; EL(:)'];  % Shape: [2 × 8] → 8 beam directions

        for i = 1:numVirtualPorts
            pos = elementPos(:, subarrayMap{i});  % [3 × subarraySize]
            a = steervec(pos, steeringDir);       % [subarraySize × 1]
            steeringVecs(subarrayMap{i},i) = mean(a,2);  % Fill entire column
        end
        
        % Normalize
        WtxVirtual = steeringVecs / norm(steeringVecs, 'fro');
        
        % Store results
        simParameters.WtxVirtual = WtxVirtual;
        simParameters.ElementPosition = elementPos;
        simParameters.SubarrayMap = subarrayMap;

        % Rx antenna array configuration in CDL channel
        rxArray = simParameters.ReceiveAntennaArray;
        M = rxArray.PanelDimensions(2);
        N = rxArray.PanelDimensions(1);
        Ng = rxArray.NumPanels;

        % channel.ReceiveAntennaArray.Size = [M N numRxPol 1 Ng];
        channel.ReceiveAntennaArray.Size = simParameters.ReceiveAntennaArray.Size;
        channel.ReceiveAntennaArray.ElementSpacing = simParameters.ReceiveAntennaArray.ElementSpacing;
        channel.ReceiveAntennaArray.PolarizationAngles = simParameters.ReceiveAntennaArray.PolarizationAngles;
        channel.ReceiveAntennaArray.Orientation = simParameters.ReceiveAntennaArray.Orientation;
        channel.ReceiveAntennaArray.Element = simParameters.ReceiveAntennaArray.Element;
        channel.ReceiveAntennaArray.PolarizationModel = simParameters.ReceiveAntennaArray.PolarizationModel;
       

    elseif contains(simParameters.DelayProfile,'TDL')

        channel = nrTDLChannel;
        channel.NumTransmitAntennas = nTxAnts;
        channel.NumReceiveAntennas = nRxAnts;

    else

        error('Channel not supported.')

    end

    % Configure common channel parameters: delay profile, delay spread, and
    % maximum Doppler shift
    channel.DelayProfile = simParameters.DelayProfile;
    channel.DelaySpread = simParameters.DelaySpread;
    channel.MaximumDopplerShift = simParameters.MaximumDopplerShift;

    % Get information about the baseband waveform after OFDM modulation step
    waveInfo = nrOFDMInfo(simParameters.Carrier);

    % Update channel sample rate based on carrier information
    channel.SampleRate = waveInfo.SampleRate;
    
end

function [channel,maxChannelDelay,N0] = setupChannel(simParameters,snrIdx)
% Reset the propagation channel. Obtain the maximum channel
% delay and calculate the AWGN standard deviation.

    % Extract channel
    channel = simParameters.Channel;
    channel.reset();

    % Get the channel information
    chInfo = info(channel);
    maxChannelDelay = chInfo.MaximumChannelDelay;

    % Calculate noise standard deviation. Normalize noise power by the IFFT
    % size used in OFDM modulation, as the OFDM modulator applies this
    % normalization to the transmitted waveform.
    SNRdB = simParameters.SNRIn(snrIdx);
    SNR = 10^(SNRdB/10);
    waveInfo = nrOFDMInfo(simParameters.Carrier);
    N0 = 1/sqrt(2.0*double(waveInfo.Nfft)*SNR);
    
    % Also normalize by the number of receive antennas if the channel
    % applies this normalization to the output
    if channel.NormalizeChannelOutputs
        N0 = N0/sqrt(chInfo.NumOutputSignals);
    end

end

function [csiReport, Hest] = initialCSIReport(simParameters,snrIdx,carrier,csirs,channel,csiFeedbackOpts)
% Get an initial CSI report using a perfect channel estimate

    [Hest,nVar] = getInitialChannelEstimate(carrier,channel,simParameters.SNRIn(snrIdx));
    Hest = virtualizeChannel(Hest, simParameters.WtxVirtual);
    csiFeedbackOpts.PerfectChannelEstimator = true;
    csirs.CSIRSPeriod = 'on';
    csiReport = hCSIEncode(carrier,csirs,Hest,nVar,csiFeedbackOpts);
    
end

function [estChannelGrid,nVar] = getInitialChannelEstimate(carrier,propchannel,SNRdB)
% Obtain channel estimate before first transmission. This can be used to
% obtain initial transmission parameters

    ofdmInfo = nrOFDMInfo(carrier);
    
    chInfo = info(propchannel);

    % Clone channel and get path gains and sample times for perfect timing
    % and channel estimation
    channel = clone(propchannel);
    release(channel);
    channel.ChannelFiltering = false;
    channel.NumTimeSamples = (ofdmInfo.SampleRate/1000/carrier.SlotsPerSubframe)+chInfo.MaximumChannelDelay;
    [pathGains,sampleTimes] = channel();
    
    % Perfect timing synch    
    pathFilters = getPathFilters(channel);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);

    % Noise variance (does not account for channel effects at this point.)
    nVar = 10^(-SNRdB/10)/size(estChannelGrid,3);
    
end

function XOverhead = getXOverhead(carrier,csirs)
% Calculate XOverhead for transport block size determination based on
% CSI-RS resource grid occupancy

    [~,csirsInfo] = nrCSIRSIndices(carrier,csirs);
    csirsRE = length(csirsInfo.KBarLBar{1})*length(csirsInfo.KPrime{1})*length(csirsInfo.LPrime{1});
    [~,XOverhead] = quantiz(csirsRE,[0 6 12],[0 6 12 18]);

end

function cdmLengths = getCSIRSCDMLengths(csirs)
%   CDMLENGTHS = getCSIRSCDMLengths(CSIRS) returns the CDM lengths given
%   the CSI-RS configuration object CSIRS.

    CDMType = csirs.CDMType;
    if ~iscell(csirs.CDMType)
        CDMType = {csirs.CDMType};
    end
    CDMTypeOpts = {'noCDM','fd-CDM2','CDM4','CDM8'};
    CDMLengthOpts = {[1 1],[2 1],[2 2],[2 4]};
    cdmLengths = CDMLengthOpts{strcmpi(CDMTypeOpts,CDMType{1})};

end

function csiFeedbackOpts = getCSIFeedbackOptions(simParameters,snrIdx)
% Create a CSI feedback algorithmic options structure

    csiFeedbackOpts = struct();
    csiFeedbackOpts.CSIReportMode = simParameters.CSIReportMode;
    csiFeedbackOpts.CSIReportPeriod = simParameters.CSIReportConfig.Period;
    csiFeedbackOpts.CSIReportConfig = simParameters.CSIReportConfig;
    csiFeedbackOpts.PerfectChannelEstimator = simParameters.PerfectChannelEstimator;
    csiFeedbackOpts.DMRSConfig = simParameters.PDSCH.DMRS;
    
    if simParameters.CSIReportMode == "AI CSI compression"
        % Copy additional link adaptation configuration for AI CSI compression mode
        csiFeedbackOpts.AINetworkFilename = simParameters.AINetworkFilename;
    
        % Download and extract a pretrained CSI network for AI CSI compression mode
        displayProgress = (snrIdx==1);
        helperCSINetDownloadData(displayProgress);
    end

end

function panelDimensions = getCSIReportPanelDimensions(antennaArray,codebookType)
% Configure the antenna array dimensions according to TS 38.214 Section
% 5.2.2.2 as a vector [N1,N2] for single-panel arrays and a vector
% [Ng,N1,N2] for multi-panel arrays.

    panelDimensions = antennaArray.PanelDimensions;
    
    % Add number of panels if codebook type is multi-panel
    if strcmpi(codebookType,'Type1MultiPanel')
        panelDimensions = [antennaArray.NumPanels panelDimensions];
    end
end

function csislot = nextCSISlot(period,nslot)
% Return the slot number of the first slot where CSI feedback can be
% reported according to the CSI report periodicity

    p = period(1); % Slot periodicity
    o = period(2); % Slot offset

    csislot = p*ceil((nslot-o)/p)+o;

end

function numElemenets = numAntennaElements(antArray)
% Calculate number of antenna elements in an antenna array

    numElemenets = antArray.NumPolarizations*antArray.NumPanels*prod(antArray.PanelDimensions);
    
end

function RIRestriction = updateRankRestriction(dmrsConfig,CSIReportConfig)
% Adjust the csi report configuration to restrict ranks, if the allowed number of DMRS Antenna ports less than the possible ranks for the given DMRS configuration  
    RIRestriction    =  CSIReportConfig.RIRestriction;
    
    if ~dmrsConfig.DMRSEnhancedR18 && (dmrsConfig.DMRSLength == 1) && strcmpi(CSIReportConfig.CodebookType, 'Type1SinglePanel')
        if (dmrsConfig.DMRSConfigurationType == 1)
            % restrict ranks upto 4 for DMRS configuration type 1
            dmrsAllowedRanks    = [ones(1,4) zeros(1,4)];
        else
            % restrict ranks upto 6 for DMRS configuration type 2
            dmrsAllowedRanks    = [ones(1,6) zeros(1,2)];
        end
        
        if isempty(RIRestriction)
            RIRestriction = ones(1,8);
        end

        RIRestriction = RIRestriction.*dmrsAllowedRanks;
    end
end

function validateCSIRSConfig(carrier,csirs,nTxAnts)
% validateCSIRSConfig(CARRIER,CSIRS,NTXANTS) validates the CSI-RS
% configuration, given the carrier specific configuration object CARRIER,
% CSI-RS configuration object CSIRS, and the number of transmit antennas
% NTXANTS.

    % Validate the number of CSI-RS ports
    if ~isscalar(unique(csirs.NumCSIRSPorts))
        error('nr5g:InvalidCSIRSPorts',...
            'All the CSI-RS resources must be configured to have the same number of CSI-RS ports.');
    end

    % Validate the CSI-RS and TX antenna array configuration
    if any(csirs.Ports_Options(csirs.RowNumber) ~= nTxAnts)
        rn = num2str(find(csirs.Ports_Options == nTxAnts),'%3d,');
        str = 'The number of CSI-RS ports must be equal to the number of Tx antenna elements. ';
        str = [str sprintf('For the Tx antenna array size configured, valid CSI-RS row numbers are (%s).',rn(1:end-1))];
        error(str)
    end

    % Validate the CDM lengths
    if ~iscell(csirs.CDMType)
        cdmType = {csirs.CDMType};
    else
        cdmType = csirs.CDMType;
    end
    if (~all(strcmpi(cdmType,cdmType{1})))
        error('nr5g:InvalidCSIRSCDMTypes',...
            'All the CSI-RS resources must be configured to have the same CDM lengths.');
    end
    if nTxAnts ~= csirs.NumCSIRSPorts(1)
        error('nr5g:InvalidNumTxAnts',['Number of transmit antennas (' num2str(nTxAnts)...
            ') must be equal to the number of CSI-RS ports (' num2str(csirs.NumCSIRSPorts(1)) ').']);
    end

    % Check for the overlap between the CSI-RS indices
    csirsInd = nrCSIRSIndices(carrier,csirs,"OutputResourceFormat",'cell');
    numRes = numel(csirsInd);
    csirsIndAll = cell(1,numRes);
    ratioVal = csirs.NumCSIRSPorts(1)/prod(getCSIRSCDMLengths(csirs));
    for resIdx = 1:numRes
        if ~isempty(csirsInd{resIdx})
            grid = nrResourceGrid(carrier,csirs.NumCSIRSPorts(1));
            [~,tempInd] = nrExtractResources(csirsInd{resIdx},grid);
            if numel(tempInd)/numel(csirsInd{resIdx}) ~= ratioVal
                error('nr5g:OverlappedCSIRSREsSingleResource',['CSI-RS indices of resource '...
                    num2str(resIdx) ' must be unique. Try changing the symbol or subcarrier locations.']);
            end
            csirsIndAll{resIdx} = tempInd(:);
            for idx = 1:resIdx-1
                overlappedInd = ismember(csirsIndAll{idx},csirsIndAll{resIdx});
                if any(overlappedInd)
                    error('nr5g:OverlappedCSIRSREsMultipleResources',['The resource elements of the '...
                        'configured CSI-RS resources must not overlap. Try changing the symbol or '...
                        'subcarrier locations of CSI-RS resource ' num2str(idx) ' and resource ' num2str(resIdx) '.']);
                end
            end
        end
    end
end

function displaySNRPointProgress(simParameters,snrIdx)
% Print SNR point progress

    str = ['\nSimulating transmission scheme 1 (%dx%d) and '...  
          'SCS=%dkHz with %s channel at %gdB SNR for %d 10ms frame(s)\n' ...
          'Using %s as CSI feedback.\n'];

    switch simParameters.CSIReportMode
        case 'RI-PMI-CQI'
            modeText = 'RI, PMI, and CQI';
        case 'AI CSI compression'
            modeText = 'compressed channel estimates';
        otherwise 
            modeText = 'channel estimates';
    end

    SNRdB = simParameters.SNRIn(snrIdx);
    fprintf(str,simParameters.NTxAnts,simParameters.NRxAnts,simParameters.Carrier.SubcarrierSpacing, ...
        simParameters.DelayProfile,SNRdB,simParameters.NFrames,modeText);

end

function printSlotInfo(NSlots,carrier,pdsch,pdschextra,blkerr,ECR,csirsTransmission,csiReports,reportIndex)
% Print information about the current slot transmission

    ncw = pdsch.NumCodewords;
    cwLayers = floor((pdsch.NumLayers + (0:ncw-1)) / ncw);
    infoStr = [];
    for cwIdx = 1:ncw
        if blkerr
            infoStrCW = "Transmission failed";
        else
            infoStrCW = "Transmission succeeded";
        end
        infoStrCW = sprintf("%22s (Layers=%d, Mod=%5s, TCR=%.3f, CR=%.3f).",infoStrCW,cwLayers(cwIdx),pdsch.Modulation,pdschextra.TargetCodeRate(cwIdx),ECR(cwIdx));
        if (ncw>1)
            infoStr = sprintf('%s\n%s%s',infoStr,sprintf('CW%d: %s',cwIdx-1),infoStrCW);
        else
            infoStr = infoStrCW;
        end
    end
    
    csirsInfoStr = [];
    if csirsTransmission
        csirsInfoStr = "CSI-RS transmission. ";
    end

    csifbInfoStr = [];
    if carrier.NSlot == 0
        csifbInfoStr = 'Using initial CSI.';
    elseif reportIndex > 0
        csifbInfoStr = sprintf("Using CSI from NSlot=%2d.",csiReports(reportIndex).NSlot);
    end

    nslot = carrier.NSlot;
    fprintf('\n(%5.2f%%) NSlot=%2d: %s',100*(nslot+1)/NSlots,nslot,join([infoStr,csifbInfoStr,csirsInfoStr]));

end


function plotLayerEVM(NSlots,nslot,pdsch,siz,pdschIndices,pdschSymbols,pdschEqSymbols)
% Plot EVM information

    persistent slotEVM;
    persistent rbEVM
    persistent evmPerSlot;
    persistent numLayers;
    
    maxNumLayers = 8;
    if (nslot==0)
        slotEVM = comm.EVM;
        rbEVM = comm.EVM;
        evmPerSlot = NaN(NSlots,maxNumLayers);
        numLayers = pdsch.NumLayers;
        figure;
    else
        % Keep the maximum number of layers in the legend
        numLayers = max(numLayers,pdsch.NumLayers);
    end
    [Ns,P] = size(pdschEqSymbols);
    pdschEqSym = zeros(Ns,maxNumLayers);
    pdschSym = zeros(Ns,maxNumLayers);
    pdschEqSym(:,1:P) = pdschEqSymbols;    
    pdschSym(:,1:P) = pdschSymbols;
    evmPerSlot(nslot+1,:) = slotEVM(pdschSym,pdschEqSym);
    subplot(2,1,1);
    plot(0:(NSlots-1),evmPerSlot,'o-');
    xlabel('Slot number');
    ylabel('EVM (%)');
    legend("layer " + (1:numLayers),'Location','EastOutside');
    title('EVM per layer per slot');

    subplot(2,1,2);
    [k,~,p] = ind2sub(siz,pdschIndices);
    rbsubs = floor((k-1) / 12);
    NRB = siz(1) / 12;
    evmPerRB = NaN(NRB,maxNumLayers);
    for nu = 1:pdsch.NumLayers
        for rb = unique(rbsubs).'
            this = (rbsubs==rb & p==nu);
            evmPerRB(rb+1,nu) = rbEVM(pdschSym(this),pdschEqSym(this));
        end
    end

    plot(0:(NRB-1),evmPerRB,'x-');
    xlabel('Resource block');
    ylabel('EVM (%)');
    legend("layer " + (1:numLayers),'Location','EastOutside');
    title(['EVM per layer per resource block, slot #' num2str(nslot)]);
    
    drawnow;
    
end

function plotCQI(simParameters,CSIReport,perc)
% Plot CQI median and percentiles 

    % Calculate median and percentiles
    med = cellfun(@(x) median([x.CQI]), CSIReport);
    p1 = cellfun(@(x) prctile([x.CQI],50-perc/2), CSIReport);
    p2 = cellfun(@(x) prctile([x.CQI],50+perc/2), CSIReport);
    
    % Calculate the percentage of CQI not in the set {median CQI -1, median CQI, median CQI +1} 
    cqiPerc = cellfun(@(x,y) sum(abs([x.CQI]-y)>1)/length(x),CSIReport,num2cell(med));
    
    figure;
    subplot(211)
    errorbar(simParameters.SNRIn,med,p2-p1,'o-.')
    ylabel('CQI value')
    title(sprintf('Median CQI and (%g,%g) Percentiles',50-perc/2,50+perc/2));
    grid

    subplot(212)
    plot(simParameters.SNRIn,cqiPerc,'o-.')
    xlabel('SNR (dB)')
    ylabel('%')
    title('CQI not in set \{median CQI -1, median CQI, median CQI +1\}');
    grid   

end