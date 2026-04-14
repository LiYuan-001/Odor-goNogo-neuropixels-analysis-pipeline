% Read all channels LFP out with specified sample ind
function dataArray = ReadBinSection(fileName,meta,Ind)

    nChan = str2double(meta.nSavedChans);

%     nFileSamp = (Ind(end)-Ind(1))*(2 * nChan);
%     
%     if round(nFileSamp)~= nFileSamp
%         error('Make sure full file is copied')
%     end
    
    %
    fid = fopen(fileName, 'rb');
    fseek(fid, (2*nChan*(Ind(1)-1)), 'bof');
    dataArray = fread(fid, [nChan,(Ind(end)-Ind(1))], 'int16=>double');
    fclose(fid);
end 