function Go_nogo_Umap(inFile,AnalyzeSes)
%
% generate 2 return + base, 2 stem + choice linearized ratemaps for each
% cell and use it for population vector decoding
% ----------------------------------------------------------------
%
% Li Yuan, Mar-25-2025

close all

addpath(genpath('C:\Users\nrpix\Documents\Matlab\umapAndEppFileExchange_4_5'));
p.timeBin = 100/10^3;
p.savePlot = 0;
p.writeToFile = 1;
p.saveDir = 'N:\yuan\AnalysisResults'; % change it to your own folder

% do umap for the run part maze and label with different way
% plot umap for this sequence
min_dist = 0.2;
spread = 2.0;
verbose = 'none';
method = 'mex';
metric = 'cosine';
see_training = false;
%                 metric = 'euclidean';
n_neighbors = 20;
n_components = 3;
% -------------------------------------------------------------------------

% Read in input information
sessInfo = SessInfoImport_Npx(inFile);

AB_color = [0,0,1;1,0,1;1,0,0;0,1,1]; % blue, magenta, red, cyan

for sessInd = AnalyzeSes(1:end)
    
    savedir = sprintf('%s%s%s%s%s',p.saveDir,'\ProcessedData\',sessInfo(sessInd).animalID,'-date-',sessInfo(sessInd).recDate);

    % load spike mat file
    spkMatName = sprintf('%s%d%s','spkMat-',p.timeBin*10^3,'ms.mat');
    spkMatFile = fullfile(savedir,spkMatName);
    spkMat = load(spkMatFile);
    % % load cell int pyr label in MEC
    % cellTypeFile = fullfile(sessInfo(i).NIDQ,'Processed', 'mec_clusterType.mat');
    % load(cellTypeFile);
    
    % clusterNum = mec_clusterType.clusterNum;
    % pyrLabel = mec_clusterType.labelInd == 1;
    % %     pyrLabel = mec_clusterType.duration >= 0.4;
    % clusterNum_Pyr = sum(pyrLabel);
    
    umap_go_nogo.animalID = sessInfo(sessInd).animalID;
    umap_go_nogo.recDate = sessInfo(sessInd).recDate;
    umap_go_nogo.odorPair1 = sessInfo(sessInd).pair1;
    umap_go_nogo.odorPair2 = sessInfo(sessInd).pair2;
    umap_go_nogo.timeBin = p.timeBin;
    
    % % imec0
    % if p.savePlot
    %     % directory for plot figures
    %     % generate a folder for each rat eah day under the current folder
    %     savedir_1 = sprintf('%s%s%d%s%s%s',cd,'\Figures\',sessInfo(i).ratID,'-day',sessInfo(i).Date,'\Umap');
    %     if ~exist(savedir_1, 'dir')
    %         mkdir(savedir_1);
    %     end
    %     delete(strcat(savedir_1,'\*'));
    % end
    % 
    

   
            spkTrain = spkMat.spkTrainMatrix.odorPair1_timeBin_spkTrain_Raw;
            spkTrain_gauss = spkMat.spkTrainMatrix.odorPair1_timeBin_spkTrain_Gauss;            
            timeBinTemp = spkMat.spkTrainMatrix.odorPair1_timeBin;
            odorLabel = spkMat.spkTrainMatrix.odorPair1_timeBin_OdorLabel;
            successLabel = spkMat.spkTrainMatrix.odorPair1_timeBin_SuccessLabel;
            outcomeLabel = spkMat.spkTrainMatrix.odorPair1_timeBin_OutcomeLabel;

            h = figure(1);
            h.Position = [100,100,900,900];    
            
            % all maze
            dataSet = spkTrain';
            %                 dataSet_z = dataSet;
            dataSet_z = zscore(dataSet,0,1);
            [reduction, umap, clusterIds, extras] = run_umap(dataSet_z,...
                'min_dist',min_dist,'metric',metric,'spread',spread,'n_components',n_components,'n_neighbors',n_neighbors, 'verbose', verbose);
    
            reduction_a = reduction(outcomeLabel ==1,:);
            reduction_b = reduction(outcomeLabel ==2,:);
            
            scatter3(reduction_a(:,1), reduction_a(:,2), reduction_a(:,3), 3, AB_color(1,:), 'filled'); % 50 = marker size, 'filled' = solid circles
            hold on
            scatter3(reduction_b(:,1), reduction_b(:,2), reduction_b(:,3), 3, AB_color(2,:), 'filled');
            % TITLE1 = 'all maze umap';
            % TITLE2 = sprintf('%s%1.2f%s%1.1f%s%s%s%d','min_dist:',min_dist,' spread:',spread,' metric:',metric,' n_neighbors:',n_neighbors);
            % title({TITLE1;TITLE2},'Interpreter','None');
            

            %% treadmill on
            subplot(2,6,2)
%             dataSet = spkMat_run(pyrLabel,delayLabel_10s_Ind & delayLabel_on_Ind)';
%             %                 dataSet_z = dataSet;
%             dataSet_z = zscore(dataSet,0,1);
%             [reduction, umap, clusterIds, extras] = run_umap(dataSet_z,...
%                 'min_dist',min_dist,'metric',metric,'spread',spread,'n_components',n_components,'n_neighbors',n_neighbors, 'verbose', verbose);
            
            valIdx = find((delayLabel_10s_Ind & delayLabel_on_Ind)==1);
            RGB_run = zeros(length(valIdx),3);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
            colormap(jet)
            % colorbar
            clim([0 30])
%             scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
            TITLE1 = 'delay on 10 s';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None');
            
            % whole 30 s
            subplot(2,6,3)
            delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
            valIdx = find((delay_30_Ind & delayLabel_on_Ind)==1);
            RGB_run = zeros(length(valIdx),3);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end          
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
            colormap(jet)
            % colorbar
            clim([0 30])
            TITLE1 = 'delay on 0 - 30 s';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None');
            
            
            % plot 10-s delay in L vs R return
            subplot(2,6,4)             
            valIdx = find((delayLabel_10s_Ind & delayLabel_on_Ind)==1);
            RGB_run = zeros(length(valIdx),3);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
%             colormap(jet)
            % colorbar
            clim([0 30])
%             scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
            TITLE1 = 'delay on 10 s ECLR';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None'); 
            
            % plot 30-s delay in L vs R return
            subplot(2,6,5)             
            delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
            valIdx = find((delay_30_Ind & delayLabel_on_Ind)==1);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
%             colormap(jet)
            % colorbar
            clim([0 30])
%             scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
            TITLE1 = 'delay on 30 s ECLR';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None'); 
            
            %% treadmill off
            subplot(2,6,8)
%             dataSet = spkMat_run(pyrLabel,delayLabel_10s_Ind & delayLabel_off_Ind)';
%             %                 dataSet_z = dataSet;
%             dataSet_z = zscore(dataSet,0,1);
%             [reduction, umap, clusterIds, extras] = run_umap(dataSet_z,...
%                 'min_dist',min_dist,'metric',metric,'spread',spread,'n_components',n_components,'n_neighbors',n_neighbors, 'verbose', verbose);
            
            valIdx = find((delayLabel_10s_Ind & delayLabel_off_Ind)==1);
            RGB_run = zeros(length(valIdx),3);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
            colormap(jet)
            % colorbar
            clim([0 30])
%             scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
            TITLE1 = 'delay off 10 s';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None');
                   

            % whole 30 s
            subplot(2,6,9)
            
            delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
            valIdx = find((delay_30_Ind & delayLabel_off_Ind)==1);
            RGB_run = zeros(length(valIdx),3);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end          
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 2, colorVal, 'filled'); % 50 = marker size, 'filled' = solid circles
            colormap(jet)
            % colorbar
            clim([0 30])
            TITLE1 = 'delay off 0 - 30 s';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None');
            
            % plot 10-s delay in L vs R return
            subplot(2,6,10)             
            valIdx = find((delayLabel_10s_Ind & delayLabel_off_Ind)==1);
            RGB_run = zeros(length(valIdx),3);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
%             colormap(jet)
            % colorbar
            clim([0 30])
%             scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
            TITLE1 = 'delay off 10 s ECLR';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None'); 
            
            % plot 30-s delay in L vs R return
            subplot(2,6,11)             
            delay_30_Ind = (delayLabel_block1_Ind + delayLabel_block2_Ind + delayLabel_block3_Ind) > 0;
            valIdx = find((delay_30_Ind & delayLabel_off_Ind)==1);
            for m = 1:length(valIdx)
                RGB_run(m,:) = ECLR_color(ECLR(valIdx(m)),:);
            end
            colorVal = delayTimeAll(valIdx);
            scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
%             colormap(jet)
            % colorbar
            clim([0 30])
%             scatter3(reduction(valIdx,1), reduction(valIdx,2), reduction(valIdx,3), 3, RGB_run, 'filled'); % 50 = marker size, 'filled' = solid circles
            TITLE1 = 'delay off 30 s ECLR';
            TITLE2 = sprintf('pyr n =: %d',sum(pyrLabel));
            title({TITLE1;TITLE2},'Interpreter','None'); 
            
            

    if p.savePlot == 1
        figure(1)
        figName = sprintf('%s%s%d%s%s%s',cd,'\Imec1-Rat-',sessInfo(sessInd).ratID,'-Day-',sessInfo(sessInd).Date,'-Fig8Run_UMAP_delay_block_time_ECLR');
        savefig(figName)
        print(figName,'-dpng','-r150');
        print(figName,'-depsc', '-painters');
        %         figName = sprintf('%s%s%d%s%s%s',savedir_1,'\Imec1-Rat-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-Fig8Run_UMAP.fig');
        %         savefig(figName)
        %         print(figName,'-dpng','-r150');
    end
    
    
    
    
    %         if p.writeToFile
    %             fileName = sprintf('%s%s%s%d%s',sessInfo(i).NIDQ,'\Processed\','rateMap_PV-',ceil(p.binWidth*10),'mm.mat');
    %             save(fileName, 'rateMap_PV');
    %         end
    close all
    fprintf('Finished ratemap analysis for session %d\n',sessInd);

end
