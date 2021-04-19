function [dicomHeader, listoffiles] = GetDicomHeader(Folder)
    % Get list of files
    current_folder = pwd;
    cd(Folder)

    listoffiles = dir('*.ima');
    numfiles = numel(listoffiles);
    if numfiles == 0
        listoffiles = dir('*.dcm');
        numfiles = numel(listoffiles);
    end

    dicomHeader = cell(1,numfiles);

    for k=1:numfiles
        filect      = listoffiles(k).name;
        dicom_info  = dicominfo(filect, 'UseDictionaryVR', true); %to avoid warnings. human_db anomysation overwrites some dicominfo parameters

        dicomHeader{k} = dicom_info;
    end
    clearvars k filect

    cd(current_folder)
end