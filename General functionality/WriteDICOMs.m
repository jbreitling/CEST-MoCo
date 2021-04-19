function WriteDICOMs(M,Slices,dicomHeader,path)
    if exist(path,'dir')
        rmdir(path,'s')
    end
    
    try
        mkdir(path);
    catch
        pause(5)
        mkdir(path);
    end

    for ii_file = 1:numel(dicomHeader)
        if ismember(dicomHeader{ii_file}.InstanceNumber, Slices) 
            M2DCM_uint16 = uint16(M(:,:,dicomHeader{ii_file}.InstanceNumber)); % select right slice
            dicomwrite(M2DCM_uint16,[path 'Slice_' num2str(dicomHeader{ii_file}.InstanceNumber) '.dcm'], dicomHeader{ii_file}, ...
                                'CompressionMode', 'none', 'CreateMode', 'copy', 'WritePrivate', true);
        end
    end
end