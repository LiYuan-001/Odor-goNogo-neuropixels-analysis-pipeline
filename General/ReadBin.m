function dataArray = ReadBin(fileName,meta,samp0, nSamp)

    nChan = str2double(meta.nSavedChans);

    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    if nargin == 2
        samp0 = 0;
        nSamp = Inf;
    else
        samp0 = max(samp0, 0);
        nSamp = min(nSamp, nFileSamp - samp0);
    end

    nSamp = min(nSamp, nFileSamp - samp0);

    sizeA = [nChan, nSamp];

    fid = fopen(fileName, 'rb');
    fseek(fid, samp0 * 2 * nChan, 'bof');
    dataArray = fread(fid, sizeA, 'int16=>double');
    fclose(fid);
end % ReadBin