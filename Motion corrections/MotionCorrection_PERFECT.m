function [M_CoReg] = MotionCorrection_PERFECT(M_Corrupted,Slices,path_target_DICOM,path_original_DICOM,T_hat)
%MotionCorrection_PERFECT  Perfect motion correction
%   Perfect motion correction using the MITK software. Uses the ground
%   truth parameters T_hat to map the motion corrupted images M_Corrupted.
%   Otherwise identical to MotionCorrection_CONV.m
%   Returns the co-registered data M_CoReg.

%% Initialization
    % Get DICOM header of moving image
    [dicomHeader, ~] = GetDicomHeader(path_original_DICOM);
    
    % Calculate transformation from the MITK reference system (scanner) to
    % the one defined with the origin being in the middle of the FOV and 
    % the axes aligned with the image axes
    T_conversion = CalculateConversionTransformation(dicomHeader);
    
    M_CoReg = NaN(size(M_Corrupted)); %initialize
    
    % Set paths for saving of temporary DICOMs
    filePath_cell = strsplit(fileparts(mfilename('fullpath')),filesep);
    filepath = strjoin(filePath_cell(1:end-1),filesep);
    path_moving_DICOM = [filepath '\moving_DICOM\'];
    path_MoCo_DICOM = [filepath '\MoCo_DICOM\']; 
    if exist(path_MoCo_DICOM,'dir')
        rmdir(path_MoCo_DICOM,'s')
    end
    mkdir(path_MoCo_DICOM);
    
    % Set algorhitm: mdra-0-13_ITKRigidSlabbedHead = Mattes mutual information
    %                XXX = XXX
    path_RegistrationAlgorhithm = [filepath '\MITK\mdra-0-13_ITKRigidSlabbedHead.dll'];

    
%% Co-Registration
% Save for each offset the DICOM, perform matching with target and read
% transformation matrix. Before the mapping the identification and
% mitigation of direct water saturation artefacts is performed.

    for ii_offset = 1
        WriteDICOMs(M_Corrupted(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
        
        % Matching
        [outcome_code, outcome_log] = system(['"' filepath '\MITK\matchR.exe'  '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm'  '" "' path_target_DICOM '" "' path_RegistrationAlgorhithm '" -o "' path_MoCo_DICOM 'RegFile.mapr']);
        if outcome_code ~= 0
            disp(outcome_log)
            error('Error occured during image registration. See log file above for more information.')
        end   
    end

    % Apply corrected transformation matrices
    for ii_offset = 1:size(M_Corrupted,4)
        WriteDICOMs(M_Corrupted(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
        
        % Save corrected transformation
        Write_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'],T_conversion*T_hat(:,:,ii_offset)*affine3d(T_conversion).invert.T)
        
        % Mapping
        [~, ~] = system(['"'  filepath '\MITK\mapR.exe' '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_target_DICOM '" -o "' path_MoCo_DICOM 'RegImage.nrrd']);
        M_CoReg(:,:,Slices,ii_offset) = nrrdread([path_MoCo_DICOM 'RegImage.nrrd']);
    end
   
%% Clean up and remove temporary directories
    rmdir(path_moving_DICOM,'s')
    rmdir(path_MoCo_DICOM,'s')
    
end