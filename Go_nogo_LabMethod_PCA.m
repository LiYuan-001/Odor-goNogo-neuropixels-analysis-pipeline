% used lab method to generate matrix and calculate PCA
% did not actively make changes for the exisiting lab code
% Li YUAN, Tohoku, 2026-May
function Go_nogo_LabMethod_PCA(inFile,analyzeSes)


p.isi_vio_Thres = 2; % percent
p.spkWidth_Thres = 0.4*10^-3; % second
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

%%  gngPCA5
%processes PCA directly from Sum file, without using matrix1 and matrix2 codes

%
%

% modified from gngSessionSum2_PCA_V2_directlyFromSumFile
% Kei Igarashi 210525


%eg) gngPCA5('N:\soma\populations\4odorsTaskResults\spike_data\Good_and_Bad_PerformanceSim1Jason_2',  1, [-3  5], [-3  4], 'optTag', 1, 'SumFileMode', 1);
%eg) gngPCA5('N:\neuroData\spike_data_sorted\LEC_AD', 1, [-1  1.2], [-1 1.3], 'mouseID', [789]); % young AD

optTag =0; % non-optTag session
% optTag =1: optTag session

% numPC =2: Uses PC1 and PC2
% numPC =3: Uses PC1, PC2 and PC3

% BehavPerformance=0 bad behavior
BehavPerformance = 1; % good behavior

SumFileMode =0;  % defalt, spikeFile will be analyzed
% SumFileMode =1:  uses sum.mat files (old format)

shuffle = 0;  % defalt, run shuffle for PCA
% shuffle =1:  load the shuffle result from

% folderName= 'N:\neuroData\spike_data_sorted\LEC_AD\YoungAD2\Layer II_III\sorted\PN'
xRange=[-2  2];
yRange=[-2 2];
numPC= 3;

% for k = 1:length(varargin)
%     if strcmpi(varargin{k},'SumFileMode')
%         SumFileMode = varargin{k+1};
%         varargin{k+1}=[];
%         varargin{k}=[];
% 
%     elseif strcmpi(varargin{k},'mouseID')
%         mouseID = varargin{k+1};
%         varargin{k+1}=[];
%         varargin{k}=[];
% 
% 
%     elseif strcmpi(varargin{k},'optTag')
%         optTag = varargin{k+1};
%         varargin{k+1}=[];
%         varargin{k}=[];
% 
%     elseif strcmpi(varargin{k},'numPC')
%         numPC = varargin{k+1};
%         varargin{k+1}=[];
%         varargin{k}=[];
% 
% 
%     elseif strcmpi(varargin{k},'shuffle')
%         shuffle = varargin{k+1};
%         varargin{k+1}=[];
%         varargin{k}=[];
%     end
% 
% end

if ~exist('SumFileMode')
    SumFileMode = 0; % defalt, spikeFile will be analyzed
end
if ~exist('optTag')
    optTag =0; % defalt, non-optTag
end
if ~exist('numPC')
    numPC = 2; % default: 2D plot of PCA
end
if ~exist('shuffle')
    shuffle = 0; % default: run shuffle
end



tic
% read in all cells
% initiate good cell label and location
cluster_depth = [];
cluster_isiVio = [];
cluster_width = [];
for sessInd = analyzeSes(1:end)

    close all
    savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate);
    savedir_spikeFile = sprintf('%s%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate,'\LabMethod\spikeFile');
     % get the calculated parameters
    clusterTypeFile = fullfile(savedir,'npx_Cluster_type.mat');
    load(clusterTypeFile);
    cluster_depth = [cluster_depth;npx_Cluster_type.imec0.ycoord];
    cluster_isiVio = [cluster_isiVio;npx_Cluster_type.imec0.isi_vio];
    cluster_width = [cluster_width;npx_Cluster_type.imec0.width];


    close all force

    if SumFileMode ==1
        fileNameRule = '*Sum.mat';
    else
        fileNameRule = 'spikeFile*.mat';
    end

    fileNameRule = 'spikeFile*.mat';
    f=dir(fullfile(savedir_spikeFile,fileNameRule));


    %% load data from all cells

    for n=1:size(f,1)

        load(fullfile(savedir_spikeFile,f(n).name))

        allCellData.odorIndex{n} =              history.odorIndex;
        allCellData.rst_all{n} =                history.rst_all;

        if optTag ==1
            allCellData.wave_corr(n,:) =            optdata.wave_corr;
            allCellData.p_SALT(n,:) =               optdata.p_SALT;
        end
        allCellData.learning_idx(n,1) =         learning_idx;
        % allCellData.spikeW_idx(n,1) =           spikeW;
        % allCellData.meanFR_idx(n,1) =           meanFR;

        allCellData.AB_PSTH_gs_A(n,:)=          history.AB_PSTH_gs{1,3};
        allCellData.AB_PSTH_gs_B(n,:)=          history.AB_PSTH_gs{2,3};


        if ~(SumFileMode ==1)
            if contains(mouse, 'mouse')
                mouse = extractAfter(mouse, 'mouse');
            end
            allCellData.mouse{n}=          mouse;

        end
    end
end


%% choose clean wide neurons  - Li YUAN
% change when needed
single_Cluster = cluster_isiVio <= p.isi_vio_Thres;
wide_Cluster = cluster_width >= p.spkWidth_Thres;
depth = cluster_depth >= 1000 & cluster_depth <= 5000;
chosenCell = single_Cluster(:) & wide_Cluster(:) & depth(:);

%% choosing optTaged cells
if optTag ==1

    %wave corr >= 0.85 (0.85 Kvistiani et al; 0.9: Cohen et al., 2012; 0.6 Nonomura et al., 2018; )
    threshold_corr=0.85;
    index_opt1=allCellData.wave_corr(:,1)>=threshold_corr;
    index_opt2=allCellData.wave_corr(:,2)>=threshold_corr;

    %SALT p-value (0.01 Kvistiani et al)
    threshold_pSALT=0.01;
    index_opt1=index_opt1 & allCellData.p_SALT(:,1)<=threshold_pSALT;
    index_opt2=index_opt2 & allCellData.p_SALT(:,2)<=threshold_pSALT;
    index_opttag = index_opt1 + index_opt2; % it allows if the cell surpass the threshold either in first OR last opttag sessions


    chosenCell =   (chosenCell ==1) & (index_opttag>0);
end

%% choosing Regular spike cell

% if optTag == 0
%    chosenCell =   (chosenCell ==1) & (allCellData.spikeW_idx(:) >= 230) & (allCellData.meanFR_idx(:) >= 0.1); %exclude interneurons
% end

%% choose cells depending on behavior performance


% % if perfm_f>=75 && perfm_n>=75
% %     learning_idx=0;disp('Good performance')
% % elseif perfm_f>=75 && perfm_n<75% familiar good
% %     learning_idx=1;disp('Bad performance for novel pairs')
% % elseif perfm_f<75 && perfm_n>=75% novelty good
% %     learning_idx=2;disp('Bad performance for familiar pairs')
% % else
% %     learning_idx=3;disp('Bad performance')
% % end


if BehavPerformance==1
    chosenCell =   (chosenCell(:) ==1) & (allCellData.learning_idx(:) == 0);% learning_idx=0: familiar and novel good

elseif BehavPerformance==0
    chosenCell =  (chosenCell(:) ==1)  & (allCellData.learning_idx(:) == 1); % learning_idx=1: familiar good but novel bad

elseif BehavPerformance==2
    chosenCell =  (chosenCell(:) ==1)  & (allCellData.learning_idx(:) >= 0); % learning_idx=0 or 1: any performance

    % elseif BehavPerformance==2
    %     chosenCell =  (chosenCell(:) ==1)  & or(allCellData.learning_idx(:) == 0,allCellData.learning_idx(:) == 1); % learning_idx=1: familiar good but novel good and bad

end




%% choose data only from specified mice

if exist('mouseID')

    mouseIDflag = zeros(size(f,1),1);

    for n=1:size(f,1)

        for k=1:size(mouseID,2)

            if mouseID(k) ==  str2num(allCellData.mouse{n})
                mouseIDflag(n) =mouseIDflag(n) + 1;

            end

        end
    end

    chosenCell = (chosenCell(:) ==1) &  (mouseIDflag(:) ==1) ;
end


%% extract data for AB session
x_range=21:141;

matrixRange.odorA1= 1:121;
matrixRange.odorB1=122:242;
matrixRange.odorC1=243:363;
matrixRange.odorD1=364:484;

matrixRange.time_odor=21:40;%0 to 1s
matrixRange.time_delay=61:80;%2 to 3s


dataMatrix_ABsession_odorA=allCellData.AB_PSTH_gs_A(chosenCell,x_range);
dataMatrix_ABsession_odorB=allCellData.AB_PSTH_gs_B(chosenCell,x_range);


CTS_data2=[dataMatrix_ABsession_odorA dataMatrix_ABsession_odorB];
CTS_data2=CTS_data2-min(CTS_data2,[],2);
CTS_data2=CTS_data2./max(CTS_data2,[],2);
CTS_data2=CTS_data2(~isnan(CTS_data2(:,1)),:); %removes cells with NaN data


%% extract data for AB12 session


odorIndex_allCells =allCellData.odorIndex(chosenCell);  %cue ID of all trials for 239 cells
rst_all_allCells =  allCellData.rst_all(chosenCell);  %raster of all trials for 239 cells


%% PCA and Figure

mode =0;
inhStimTag = 0;
firstLast =2;

toc

savedir_PCA = sprintf('%s%s',p.saveDir,'\LabMethod');
dateStr = datestr(now, 'yyyymmdd_HHMM');
sesStr = sprintf('%d_', analyzeSes);
sesStr = sesStr(1:end-1);   % remove final "_"

folderName = fullfile(savedir_PCA, 'PCARun', sprintf('%s_ses%s', dateStr, sesStr));
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
% Save session info text file
infoFile = fullfile(folderName, 'analyzeSes_info.txt');
fid = fopen(infoFile, 'w');
fprintf(fid, 'analyzeSes:\n');
fprintf(fid, '%s\n\n', mat2str(analyzeSes));
fprintf(fid, 'Session details:\n');
fprintf(fid, 'SessionIndex\tAnimalID\tRecDate\n');
for i = 1:numel(analyzeSes)
    sesInd = analyzeSes(i);
    fprintf(fid, '%d\t%s\t%s\n', ...
        sesInd, ...
        string(sessInfo(sesInd).animalID), ...
        string(sessInfo(sesInd).recDate));
end
fclose(fid);

%%AB-only session
[meanFivePercentile_AB] = trialShufflePCA(folderName,odorIndex_allCells, rst_all_allCells, matrixRange, x_range, numPC, mode, BehavPerformance,inhStimTag,firstLast,shuffle);
[ISI_ABsession] = gngSessionSum2_PCA_AB2_v3(CTS_data2, matrixRange,numPC,xRange, yRange, mode,BehavPerformance,meanFivePercentile_AB) ;
%% ABCD session
meanFivePercentile=[];
[meanFivePercentile] = trialShufflePCA_ABCD(folderName,odorIndex_allCells, rst_all_allCells, matrixRange, x_range, numPC, mode, BehavPerformance,inhStimTag,firstLast,shuffle);
%bootstrapPCA(odorIndex_allCells, rst_all_allCells, Time,x_range, numPC, mode, BehavPerformance,inhStimTag,firstLast,meanFivePercentile);
gngPCA_T1_T5_250411(odorIndex_allCells, rst_all_allCells, matrixRange, x_range, numPC, meanFivePercentile,inhStimTag,mode,BehavPerformance, ISI_ABsession, xRange, yRange);
%

% print(figName,'-dpng','-r200');
% savefig(gcf, 'PCA-wholeHPC-91cells.fig')
end








