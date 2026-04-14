% Read single channel LFP out
function dataArray = ReadBinChannel(fileName,meta,ch)

    nChan = str2double(meta.nSavedChans);

    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    
    if round(nFileSamp)~= nFileSamp
        error('Make sure full file is copied')
    end
    
    % skip 384 channel data
    skip = 2*(nChan-1);
    fid = fopen(fileName, 'rb');
    fseek(fid, 2 * ch-2, 'bof');
    dataArray = fread(fid, nFileSamp, 'int16=>double',skip);
    fclose(fid);
end