function [M_CoReg, T] = MotionCorrection_PROPOSED(M_Corrupted,Slices,path_target_DICOM,path_original_DICOM)
%MotionCorrection_PROPOSED  Conventional motion correction with the
%additional identification and mitigation of direct water saturation
%artefacts
%   A conventional motion correction using the MITK software with the
%   additional identification and mitigation of direct water saturation
%   artefacts. Co-registers the motion corrupted images M_Corrupted
%   (x,y,z,dw) using only the slices specified to the target image stored
%   in the folder path_target_DICOM. As MITK requires DICOMs the folder of
%   an original DICOM is required: path_original_DICOM.
%   Returns the co-registered data M_CoReg and the 4x4 transformation
%   matrices T(:,;,dw). T are defined with the origin being in the middle
%   of the FOV and the axes aligned with the image axes.

%% Initialization
    % Get DICOM header of moving image
    [dicomHeader, ~] = GetDicomHeader(path_original_DICOM);
    
    % Calculate transformation from the MITK reference system (scanner) to
    % the one defined with the origin being in the middle of the FOV and 
    % the axes aligned with the image axes
    T_conversion = CalculateConversionTransformation(dicomHeader);
    
    M_CoReg = NaN(size(M_Corrupted)); %initialize
    T = NaN(4,4,size(M_Corrupted,4)); %initialize
    
    % Set paths for saving of temporary DICOMs
    filePath_cell = strsplit(fileparts(mfilename('fullpath')),filesep);
    filepath = strjoin(filePath_cell(1:end-1),filesep);
    path_moving_DICOM = [filepath '\moving_DICOM_'];
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

      for ii_offset = 1:size(M_Corrupted,4)
        
        WriteDICOMs(M_Corrupted(:,:,:,ii_offset),Slices,dicomHeader,[path_moving_DICOM num2str(ii_offset) '\']);
        
        % Matching
        [outcome_code, outcome_log] = system(['"' filepath '\MITK\matchR.exe'  '" "' path_moving_DICOM num2str(ii_offset) '\Slice_' num2str(Slices(1)) '.dcm'  '" "' path_target_DICOM '" "' path_RegistrationAlgorhithm '" -o "' path_MoCo_DICOM 'RegFile.mapr']);
        if outcome_code ~= 0
            disp(outcome_log)
            error('Error occured during image registration. See log file above for more information.')
        end
        
        % Read transformation matrix
        T(:,:,ii_offset) = affine3d(T_conversion).invert.T*Read_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'])*T_conversion;    
    end

    
    % Identification and mitigation of direct water saturation artefacts
    % Weights wk for each offset calculated as L2-norm
    wk = squeeze(sqrt(sum(M_Corrupted(:,:,Slices,:).^2,[1 2 3])));
    T = AIM(T,wk.^2);
    
    % Apply corrected transformation matrices
    for ii_offset = 1:size(M_Corrupted,4)
        %WriteDICOMs(M_Corrupted(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
        
        % Save corrected transformation
        Write_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'],T_conversion*T(:,:,ii_offset)*affine3d(T_conversion).invert.T)
        
        % Mapping
        [~, ~] = system(['"'  filepath '\MITK\mapR.exe' '" "' path_moving_DICOM num2str(ii_offset) '\Slice_' num2str(Slices(1)) '.dcm' '" "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_target_DICOM '" -o "' path_MoCo_DICOM 'RegImage.nrrd']);
        M_CoReg(:,:,Slices,ii_offset) = nrrdread([path_MoCo_DICOM 'RegImage.nrrd']);
    end
   
%% Clean up and remove temporary directories
    for ii_offset = 1:size(M_Corrupted,4)
        rmdir([path_moving_DICOM num2str(ii_offset)],'s')
    end
    rmdir(path_MoCo_DICOM,'s')
    
end

function T = AIM(T,wk)
    
%% Parameters
    % Radius specifying the volume of interest. For head empirically
    % determined as 70 mm.
    R = 70; %[mm]
    % Standard deviation of Gaussian used for weighting of neighboring
    % measurements.
    sigma = 0.6;

%% Identification and mitigation (Figure 2)
    bool_thereAreOutliers = true;
    while bool_thereAreOutliers
        d_RMS_i = NaN(size(T,3),1); %initialize
        T_tilde = repmat(zeros(4),1,1,size(T,3)); %initialize

        for ii_offset = 1:size(T,3)
            
            % Calculate weights. Note, wk is normalized in the exp function
            % by the largest wk to avoid numerical instabilities.
            wi = wk'.*exp(-(((1:size(T,3)) - ii_offset).^2./2/(sigma/(wk(ii_offset)/max(wk)+0.2))^2));
            wi = wi/sum(wi);

            % Log-matrix blending (Equation XXX)
            for kk_offset = 1:size(T,3)
                T_tilde(:,:,ii_offset) = T_tilde(:,:,ii_offset) + wi(kk_offset)*logm(T(:,:,kk_offset));
            end
            
            T_tilde(:,:,ii_offset) = expm(T_tilde(:,:,ii_offset));
            T_tilde(:,4,ii_offset) = [0 0 0 1]; %Required to avoid small numerical errors.
            
            % Root-mean-square deviation as difference measure (Equation XXX)
            M = T(:,:,ii_offset)/T_tilde(:,:,ii_offset) - eye(4);
            d_RMS_i(ii_offset) = sqrt(1/5 * R^2 * trace((M(1:3,1:3))'*M(1:3,1:3)) + (M(4,1:3))*(M(4,1:3))');
        end
        
        % Identify outliers using MAD
        bool_outlier = isoutlier(d_RMS_i,'ThresholdFactor',10);
        
        % If there is an outlier, replace largest with weighted average.
        % Else AIM is finished.
        if (sum(bool_outlier) >0)
            [~,max_outlier] = max(d_RMS_i);
            T(:,:,max_outlier) = T_tilde(:,:,max_outlier);
        else
            bool_thereAreOutliers = false;
        end
    end
end