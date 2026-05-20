inFile = 'probe_analysis_LiYUAN.xlsx';
analyzeSes = [2,6:10];

% compulsary preprocess by order
Preprocess_Allclocks_Sync(inFile,analyzeSes) % shared by whoever use this data, only need to run it once for the dataset
Preprocess_sessionSplit2(inFile,analyzeSes) % shared
Preprocess_Behavior(inFile,analyzeSes)
Preprocess_Npx_extractClusters(inFile,analyzeSes)
Go_nogo_LabMethod_matFile_genetation(inFile,analyzeSes)


% % optional cell analysis
% Go_nogo_cell2task_LabMethod(inFile,analyzeSes)
% Go_nogo_cell2taskRateMap(inFile,analyzeSes)
