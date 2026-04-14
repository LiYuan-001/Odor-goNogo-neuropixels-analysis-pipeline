% check raw waveform by the spikeGLX and decide the time used for
% visualization
% Li YUAN, 2026-Mar-09

in_binFile = "D:\m00005\m00005_po0_20260402_headfix-2_g0\m00005_po0_20260402_headfix-2_g0_imec0\m00005_po0_20260402_headfix-2_g0_t0.imec0.ap.bin";
in_binFile2 = "D:\m00005\m00005_po0_20260402_headfix-2_g0\catgt_m00005_po0_20260402_headfix-2_g0_imec0\m00005_po0_20260402_headfix-2_g0_tcat.imec0.ap.bin";
in_metaFile = "D:\m00005\m00005_po0_20260402_headfix-2_g0\m00005_po0_20260402_headfix-2_g0_imec0\m00005_po0_20260402_headfix-2_g0_t0.imec0.ap.meta";
% chanSelect = [116:-2:96,46:-2:0,334:-2:328,324:-2:288];
% chanSelect = [264:-4:240,188:-4:144,92:-4:54];
chanSelect = [285,281,277,273,269];
timeSelect = [900.6];
p.timeRange = 0.2;  % unit sec

% set parameter
p.savePlot = 1;
p.theta = [6 12];
p.gamma = [30 90];
p.yGap = 90;
MmV = 400;
MmV_gamma = 100;


% MmV = max(max(abs(dataArray_detrend))); % absolute maximum voltage

% Read in input information

close all

% 
% if p.savePlot
%     % directory for plot figures
%     % generate a folder for each rat eah day under the current folder
%     savedir_1 = sprintf('%s%s%d%s%s%s',cd,'\Figures\LFP_surveyProbe_CSDcheck');
%     if ~exist(savedir_1, 'dir')
%         mkdir(savedir_1);
%     end
% end
    
        % load file
        meta = ReadMeta(in_metaFile);
        converter = ADconvert(meta); % convert to Volt
        Fs = 30000;
        if isfield(meta, 'imDatPrb_type')
            pType = str2num(meta.imDatPrb_type);
        else
            pType = 0; %3A probe
        end
        
 
        % load channel map (this is default first bank)
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta);
        chanInd = chanSelect + 1;
        yDepth = yCoord(chanInd);


        for m = 1:length(timeSelect)
            h = figure;
            h.Position = [100,100,1200,900];
            sampleInd = ceil(timeSelect(m)*Fs);
            lfpInd = sampleInd:(sampleInd + p.timeRange * Fs);
            dataArray = ReadBinSection(in_binFile,meta,lfpInd)*converter.AP(1)*10^6; % unit uV
            dataArray2 = dataArray(chanInd,:);

            dataArray_spikeTemp = ReadBinSection(in_binFile2,meta,lfpInd)*converter.AP(1)*10^6; % unit uV
            dataArray_Spike = dataArray_spikeTemp(chanInd,:);

            


            % dataArray_Spike = zeros(size(dataArray2));
            % dwonsample to 2500;
            Fs_LFP = 2500;
            dataArray3 = resample(dataArray2', Fs_LFP, Fs);
            dataArray_detrend = detrend(dataArray3);
            dataArray_detrend = dataArray_detrend';
            lfplength = size(dataArray3,1);
            dataArray_LFP = zeros(size(dataArray3'));
            dataArray_Ripple = zeros(size(dataArray3'));

            for nn = 1:length(chanInd)
                dataArray_LFP(nn,:) = fftbandpass(dataArray_detrend(nn,:),Fs_LFP,1,2,199,200);
                dataArray_Ripple(nn,:) = fftbandpass(dataArray_detrend(nn,:),Fs_LFP,100,101,249,250);
            end

            subplot(2,1,1)
            for nn = 1:length(chanInd)
                plot((1:(lfplength))/Fs_LFP,dataArray_LFP(nn,:)-600*nn,'k','LineWidth',1)
                hold on
            end

            subplot(2,1,2)
            for nn = 1:length(chanInd)
                plot((1:length(lfpInd(2:end)))/Fs,dataArray_Spike(nn,:)-600*nn,'k','LineWidth',1)
                hold on
            end
            


            % subplot(1,3,2)
            % for nn = 1:length(chanInd)
            %     plot((1:(lfplength))/Fs_LFP,dataArray_LFP(nn,:)/MmV*15+yDepth(nn))
            %     hold on
            % end

            for nn = 1:length(chanInd)
                plot((1:(lfplength))/Fs_LFP,dataArray_LFP(nn,:)-800*nn,'k','LineWidth',2)
                text(-0.07,-800*nn,num2str(yDepth(nn)))
                hold on
                xlim([-0.08 1])
            end

            % subplot(1,3,3)
            % for nn = 1:length(chanInd)
            %     plot((1:(lfplength))/Fs_LFP,dataArray_Ripple(nn,:)/MmV*15+yDepth(nn))
            %     hold on
            % end

        end

      
            
            % initiate the variables for the phase and amp
            dataArray_Theta = zeros(size(dataArray3));
            dataArray_Gamma = zeros(size(dataArray3));
            dataArray_detrend = detrend(dataArray3);
            
            phase_Theta = zeros(size(dataArray3)); 
            amp_Theta = zeros(size(dataArray3)); 
            phase_Gamma = zeros(size(dataArray3)); 
            amp_Gamma = zeros(size(dataArray3)); 
            
            for nn = 1:length(chanInd)
                dataArray_Theta(:,nn) = fftbandpass(dataArray_detrend(:,nn),Fs_LFP,p.theta(1)-1,p.theta(1),p.theta(2),p.theta(2)+1);
                dataArray_Gamma(:,nn) = fftbandpass(dataArray_detrend(:,nn),Fs_LFP,p.gamma(1)-1,p.gamma(1),p.gamma(2),p.gamma(2)+1);
                
                % [phase_Theta(:,nn),amp_Theta(:,nn)] = thetaPhase2(dataArray_Theta(:,nn));  
                % [phase_Gamma(:,nn),amp_Gamma(:,nn)] = thetaPhase2(dataArray_Gamma(:,nn));  
%                 
%                 phase_Theta(:,nn) = rad2deg(phase_Theta(:,nn));
%                 phase_Gamma(:,nn) = rad2deg(phase_Gamma(:,nn));
            end
                    
            % % plot channel map
            % subplot(2,3,1)
            % xCoord2 = xCoord + shankInd2(m)*xshankStep;
            % yCoord2 = yCoord;           
            % plotSaved(xCoord2, yCoord2, shankInd, meta)
            % TITLE1 = animal_Session;
            % TITLE2 = sprintf('%s%d%s%d%s%d','Probe: ',probeID(k),' Shank: ',shankInd2(m)+1, ' Bank ',shankBank2(m)+1);
            % title({TITLE1;TITLE2},'Interpreter','None')
            % 
            % ycoordSelect = yCoord2(chanInd);
        
            % calculate CSD based on the selected LFP and plot CSD and LFP,
            % Theta
            f(1) = subplot(2,3,2);
            [CSDoutput]  = CSD_simple(dataArray_Theta); % unit: volt, meter
            xAxis = timeSelect(m) + (1:lfplength)./Fs_LFP;
            cmax = max(max(CSDoutput));
            contourf(xAxis,ycoordSelect(2:end-1),CSDoutput',40,'LineColor','none');hold on;
            colormap jet; 
            caxis([-cmax cmax]/2);
            set(gca,'YDir','normal');
            xlabel('time (s)');ylabel('depth');title('CSD');
            hold on
            
            for nn = 1:length(chanInd)
                plot(xAxis,dataArray_detrend(:,nn)/MmV*p.yGap+ycoordSelect(nn),'Color',[0.3,0.3,0.3])
                hold on
            end
            TITLE1 = sprintf('%s%d%s','Theta range CSD, LFP space: ',p.yGap,'um');
            title({TITLE1},'Interpreter','None')
            
            % gamma
            f(2) = subplot(2,3,3);
            plotRange = 1:round(Fs_LFP/5);
            [CSDoutput]  = CSD_simple(dataArray_Gamma); % unit: volt, meter
            xAxis = timeSelect(m) + (1:lfplength)./Fs_LFP;
            cmax = max(max(CSDoutput));
            contourf(xAxis(plotRange),ycoordSelect(2:end-1),CSDoutput(plotRange,:)',40,'LineColor','none');hold on;
            colormap jet; 
            caxis([-cmax cmax]/2);
            set(gca,'YDir','normal');
            xlabel('time (s)');ylabel('depth');title('CSD');
            hold on
            
            for nn = 1:length(chanInd)
                plot(xAxis(plotRange),dataArray_Gamma(plotRange,nn)/MmV_gamma*p.yGap+ycoordSelect(nn),'Color',[0.3,0.3,0.3])
                hold on
            end
            TITLE1 = sprintf('%s','Gamma range CSD');
            title({TITLE1},'Interpreter','None')
            
            linkaxes(f,'y')
            
            % calculate power, phase shift and coherence
            subplot(2,3,4)
            theta_pow_mean = nanmean(amp_Theta,1);
            theta_pow_err = std(amp_Theta,[],1)./sqrt(size(amp_Theta,1));
            gamma_pow_mean = nanmean(amp_Gamma,1);
            gamma_pow_err = std(amp_Gamma,[],1)./sqrt(size(amp_Gamma,1));
%             plot(ycoordSelect,theta_pow_mean,'k')
%             hold on
            errorbar(ycoordSelect,theta_pow_mean,theta_pow_err,theta_pow_err,'k')
            view(90,-90)
            hold on
            errorbar(ycoordSelect,gamma_pow_mean,gamma_pow_err,gamma_pow_err,'r')
            xlabel('channel depth (um)')
            ylabel('Power')
            title('Theta power: black;  Gamma power: red')
            
%             subplot(2,6,9)
%             theta_phaseZero = phase_Theta(:,1);
%             theta_phaseDiff = phase_Theta - theta_phaseZero;
%             theta_phase_mean = rem(nanmean(theta_phaseDiff,1),360);
%             theta_phase_err = std(theta_phaseDiff,[],1)./sqrt(size(theta_phaseDiff,1));
%             
%             gamma_phaseZero = phase_Gamma(:,1);
%             gamma_phaseDiff = phase_Gamma - gamma_phaseZero;
%             gamma_phase_mean = nanmean(gamma_phaseDiff,1);
%             gamma_phase_err = std(gamma_phaseDiff,[],1)./sqrt(size(gamma_phaseDiff,1));
%             
%             errorbar(ycoordSelect,theta_phase_mean,theta_phase_err,theta_phase_err,'k')
%             view(90,-90)
%             xlabel('channel depth (um)')
%             ylabel('Phase diff (deg)')
%             title('Theta phase shift')
%             
%             subplot(2,6,10)
%             errorbar(ycoordSelect,gamma_phase_mean,gamma_phase_err,gamma_phase_err,'r')
%             view(90,-90)
%             xlabel('channel depth (um)')
%             ylabel('Phase diff (deg)')
%             title('Gamma phase shift')
            
            
            % % save plot photots
            % if p.savePlot == 1
            %     figName = sprintf('%s%s%s%s%d%s%d%s%d%s',savedir_1,'\',animal_Session,'-probe-',probeID(k),...
            %         '-shank-',shankInd2(m)+1, '-Bank-',shankBank2(m)+1,'-LFP-CSD-shift');
            %     print(figName,'-dpng','-r300');
            % end           


    fprintf('Finished reading raw signals\n');
    close all


% =========================================================
% Parse snsGeomMap for XY coordinates
%
function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)

    C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
    shankInd = double(cell2mat(C(1)));
    xCoord = double(cell2mat(C(2)));
    yCoord = double(cell2mat(C(3)));
    connected = double(cell2mat(C(4)));

    % parse header for number of shanks
    geomStr = meta.snsGeomMap;
    headStr = extractBefore(geomStr,')(');
    headParts = split(headStr,',');
    nShank = str2double(headParts{2});
    shankWidth = str2double(headParts{4});
    shankPitch = str2double(headParts{3});
end % geomMapToGeom


function geom = getGeomParams(meta)
% create map
geomTypeMap = makeTypeMap();

% get probe part number; if absent, this is a 3A
if isfield(meta,'imDatPrb_pn')
    pn = meta.imDatPrb_pn;
else
    pn = '3A';
end

if geomTypeMap.isKey(pn)
    geomType = geomTypeMap(pn);
else
    fprintf('unsupported probe part number\n');
    return;
end

switch geomType
    case 'np1_stag_70um'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'nhp_lin_70um'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'nhp_stag_125um_med'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 87;
        geom.vertPitch = 20;
        geom.rowsPerShank = 1368;
        geom.elecPerShank = 2496;
    case 'nhp_stag_125um_long'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 87;
        geom.vertPitch = 20;
        geom.rowsPerShank = 2208;
        geom.elecPerShank = 4416;
    case 'nhp_lin_125um_med'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 103;
        geom.vertPitch = 20;
        geom.rowsPerShank = 1368;
        geom.elecPerShank = 2496;
    case 'nhp_lin_125um_long'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 103;
        geom.vertPitch = 20;
        geom.rowsPerShank = 2208;
        geom.elecPerShank = 4416;
    case 'uhd_8col_1bank'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 14;
        geom.odd_xOff = 14;
        geom.horzPitch = 6;
        geom.vertPitch = 6;
        geom.rowsPerShank = 48;
        geom.elecPerShank = 384;
   case 'uhd_8col_16bank'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 14;
        geom.odd_xOff = 14;
        geom.horzPitch = 6;
        geom.vertPitch = 6;
        geom.rowsPerShank = 768;
        geom.elecPerShank = 6144;
    case 'np2_ss'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 15;
        geom.rowsPerShank = 640;
        geom.elecPerShank = 1280;
     case 'np2_4s'
        geom.nShank = 4;
        geom.shankWidth = 70;
        geom.shankPitch = 250;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 15;
        geom.rowsPerShank = 640;
        geom.elecPerShank = 1280;
    case 'NP1120'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 6.75;
        geom.odd_xOff = 6.75;
        geom.horzPitch = 4.5;
        geom.vertPitch = 4.5;
        geom.rowsPerShank = 192;
        geom.elecPerShank = 384;
    case 'NP1121'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 6.25;
        geom.odd_xOff = 6.25;
        geom.horzPitch = 3;
        geom.vertPitch = 3;
        geom.rowsPerShank = 384;
        geom.elecPerShank = 384;
    case 'NP1122'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 12.5;
        geom.odd_xOff = 12.5;
        geom.horzPitch = 3;
        geom.vertPitch = 3;
        geom.rowsPerShank = 24;
        geom.elecPerShank = 384;
    case 'NP1123'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 10.25;
        geom.odd_xOff = 10.25;
        geom.horzPitch = 4.5;
        geom.vertPitch = 4.5;
        geom.rowsPerShank = 32;
        geom.elecPerShank = 384;
    case 'NP1300'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 48;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'NP1200'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 64;
        geom.elecPerShank = 128;
    case 'NXT3000'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 53;
        geom.odd_xOff = 53;
        geom.horzPitch = 0;
        geom.vertPitch = 15;
        geom.rowsPerShank = 128;
        geom.elecPerShank = 128;
    otherwise
        % shouldn't see this case
        fprintf('unsupported probe part number\n');
        return;
end
end % getGeomParams

% =========================================================
% Return geometry paramters for supported probe types
% Note that geom only contains enough info to calculate
% positions for the electrodes listed in snsShankMap
%
function M = makeTypeMap()
% many part numbers have the same geometry parameters ;
% make a map that pairs geometry type (value) with probe part number (key)
M = containers.Map('KeyType','char','ValueType','char');

M('3A') = 'np1_stag_70um';
M('PRB_1_4_0480_1') = 'np1_stag_70um';
M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
M('NP1010') = 'np1_stag_70um'; 
M('NP1011') = 'np1_stag_70um';
M('NP1012') = 'np1_stag_70um';
M('NP1013') = 'np1_stag_70um';

M('NP1015') = 'nhp_lin_70um';
M('NP1015') = 'nhp_lin_70um';
M('NP1016') = 'nhp_lin_70um';
M('NP1017') = 'nhp_lin_70um';
   
M('NP1020') = 'nhp_stag_125um_med';
M('NP1021') = 'nhp_stag_125um_med';
M('NP1030') = 'nhp_stag_125um_long';
M('NP1031') = 'nhp_stag_125um_long';

M('NP1022') = 'nhp_lin_125um_med';
M('NP1032') = 'nhp_lin_125um_long';

M('NP1100') = 'uhd_8col_1bank';
M('NP1110') = 'uhd_8col_16bank';

M('PRB2_1_2_0640_0') = 'np2_ss';
M('PRB2_1_4_0480_1') = 'np2_ss';
M('NP2000') = 'np2_ss';
M('NP2003') = 'np2_ss';
M('NP2004') = 'np2_ss';

M('PRB2_4_2_0640_0') = 'np2_4s';
M('PRB2_4_4_0480_1') = 'np2_4s';
M('NP2010') = 'np2_4s';
M('NP2013') = 'np2_4s';
M('NP2014') = 'np2_4s';

M('NP1120') = 'NP1120';
M('NP1121') = 'NP1121';
M('NP1122') = 'NP1122';
M('NP1123') = 'NP1123';
M('NP1300') = 'NP1300';

M('NP1200') = 'NP1200';
M('NXT3000') = 'NXT3000';
end % makeTypeMap

% =========================================================
% Plot x z positions of all electrodes and saved channels
%
function plotSaved(xCoord, yCoord, shankInd, meta)
    % get geometry
    g = getGeomParams(meta);

    % calculate positions on one shank
    nCol = g.elecPerShank/g.rowsPerShank;
    rowInd = 0:g.elecPerShank-1;
    rowInd = floor(rowInd./nCol);
    oddRows = logical(mod(rowInd,2));
    evenRows = ~oddRows;

    colInd = 0:g.elecPerShank-1;
    colInd = mod(colInd,nCol);

    xall = colInd*g.horzPitch;
    xall(evenRows) = xall(evenRows) + g.even_xOff;
    xall(oddRows) = xall(oddRows) + g.odd_xOff;

    yall = rowInd*g.vertPitch;
    
%     figure('Name','shank view','Units','Normalized', 'Position', [0.2,0.1,0.3,0.8])
    for sI = 0:g.nShank-1
        cc = find(shankInd == sI);
        scatter( g.shankPitch*sI + xall, yall, 3, 'k', 'square' ); hold on;
        scatter( g.shankPitch*sI + xCoord(cc), yCoord(cc), 30, 'green', 'square', 'filled' ); hold on; 
    end
    xlim([min(xall)-5, (g.nShank-1)*g.shankPitch + max(xall)+5]);
    ylim([min(yall)-5, max(yall)+5]);
    hold off;
end % plotSaved
