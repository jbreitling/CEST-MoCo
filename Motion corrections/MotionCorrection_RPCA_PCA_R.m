function [C, T] = MotionCorrection_RPCA_PCA_R(A,Slices,path_original_DICOM,lambda_rel,d_PCA)
%MotionCorrection_RPCA_PCA_R  Cutting-edge motion correction using a 
% two-stage registration scheme using RPCA and PCA
%   Two-stage registration scheme using RPCA and PCA as proposed in:
%
%   'Bie C, Liang Y, Zhang L, et al. Motion Correction Of Chemical Exchange
%   Saturation Transfer MRI Series Using Robust Principal Component Analysis
%   (RPCA) And PCA. Quant Imaging Med Surg. 2019;9(10):1697-1713.'
%
%   The variable nominclature follows the original paper.The algorithm 
%   co-registers the motion corrupted images A (x,y,z,dw) using only the
%   slices specified to each other. As MITK requires DICOMs the folder of
%   an original DICOM is required: path_original_DICOM.
%   Returns the co-registered data C and the 4x4 transformation
%   matrices T(:,;,dw). T are defined with the origin being in the middle
%   of the FOV and the axes aligned with the image axes. In contrast to the
%   original implementation, the algorhitm sets the first image to being
%   unmoved as the iterative nature would otherwise result in a new
%   coordinate system.

    [dicomHeader, ~] = GetDicomHeader(path_original_DICOM);
    T_conversion = CalculateConversionTransformation(dicomHeader);
    
    filePath_cell = strsplit(fileparts(mfilename('fullpath')),filesep);
    filepath = strjoin(filePath_cell(1:end-1),filesep);
    path_moving_DICOM = [filepath '\moving_DICOM\'];
    path_target_DICOM = [filepath  '\target_DICOM\'];
    path_MoCo_DICOM = [filepath  '\MoCo_DICOM\'];
    path_RegistrationAlgorhithm = [filepath  '\MITK\mdra-0-13_ITKRigidSlabbedHead.dll'];
    path_Mask = [filepath '\Mask.nrrd'];
    path_Mask_target = [filepath  '\Mask_target.nrrd'];
    
    if exist(path_MoCo_DICOM,'dir')
        rmdir(path_MoCo_DICOM,'s')
    end
    mkdir(path_MoCo_DICOM);
    
    T = repmat(eye(4),1,1,size(A,4));
    
    R_spatial = imref3d(size(A(:,:,:,1)),...
        size(A,2)*dicomHeader{1}.PixelSpacing(2).*[-0.5 +0.5],...
        size(A,1)*dicomHeader{1}.PixelSpacing(1).*[-0.5 +0.5],...
        size(A,3)*dicomHeader{1}.SliceThickness.*[-0.5 +0.5]);

    Mask_orig = zeros(size(A(:,:,:,1)));
    Mask_orig(:,:,Slices) = 1;
    Mask_offset = repmat(Mask_orig,1,1,1,size(A,4));
    
    %% STAGE I: Coarse registration using the average of RPCA
    B = A;
    
    Mtemp = permute(A,[4 1 2 3]);
    M1 = Mtemp(:,Mask_orig == 1)';
    
    lambda_0 = 1/sqrt(max(size(M1)));
    
    % RPCA algorithm taken from https://github.com/dlaptev/RobustPCA
    [L_star,~] = RobustPCA(M1,lambda_rel*lambda_0);
    
    PCi = zeros(size(Mtemp));
    PCi(:,Mask_orig == 1) = L_star';
    PCi = permute(PCi,[2 3 4 1]);
    
    average_L = mean(PCi,4);
    
    WriteDICOMs(average_L,Slices,dicomHeader,path_target_DICOM);
    
    for ii_offset = 1:size(A,4)
        WriteDICOMs(A(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
        [outcome_code, outcome_log] = system(['"' filepath '\MITK\matchR.exe'  '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm'  '" "' path_target_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_RegistrationAlgorhithm '" -o "' path_MoCo_DICOM 'RegFile.mapr' '"']);
        if outcome_code ~= 0
            disp(outcome_log)
            warning('Error occured during image registration. See log file above for more information.\n')
        end
        T(:,:,ii_offset) = affine3d(T_conversion).invert.T*Read_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'])*T_conversion;
    end

    for ii_offset = size(A,4):-1:1
        T(:,:,ii_offset) = T(:,:,ii_offset)*affine3d(T(:,:,1)).invert.T;
        WriteDICOMs(A(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
        Write_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'],T_conversion*T(:,:,ii_offset)*affine3d(T_conversion).invert.T)
        [~, ~] = system(['"' filepath  '\MITK\mapR.exe' '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_target_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" -o "' path_MoCo_DICOM 'RegImage.nrrd' '"']);
        B(:,:,Slices,ii_offset) = nrrdread([path_MoCo_DICOM 'RegImage.nrrd']);
    end
    
    %% STAGE II: Refinement registration using PCA
    
    C = B;
    
    for ii_iteration = 1:d_PCA
        
        for ii_offset =  1:size(A,4)
            Mask_offset(:,:,:,ii_offset) = imwarp(Mask_orig, R_spatial, invert(affine3d(T(:,:,ii_offset))), 'FillValues', 0, 'OutputView', R_spatial);
        end
        
        if ii_iteration == 1
            Mask_target = Mask_offset;
            T_coarse = T;
        end
        
        Mask_ii = (prod(Mask_offset,4) == 1);
        
        if sum(Mask_ii,'all') == 0
            warning(['Mask of RPCA+PCA_R algorithm is empty. No more overlap between different offsets. Algorithm was stopped after ' num2str(ii_iteration-1) ' iterations!'])
            
            % Clean up and remove temporary directories
            rmdir(path_moving_DICOM,'s')
            rmdir(path_MoCo_DICOM,'s')
            rmdir(path_target_DICOM,'s')
            if d_PCA >0
                delete(path_Mask)
                delete(path_Mask_target)
            end
            break
        end
        
        WriteMask(dicomHeader,double(Mask_ii),path_Mask);
        
        Mtemp = permute(C,[4 1 2 3]);
        M2 = Mtemp(:,Mask_ii)';
        
        [Ui,~] = eigs(cov(M2),ii_iteration);
        
        Mi= (M2*Ui)*Ui';

        PCi = zeros(size(Mtemp));
        PCi(:,Mask_ii) = Mi';
        PCi = permute(PCi,[2 3 4 1]);
        
        for ii_offset = 1:size(B,4)
            
            WriteMask(dicomHeader,double(Mask_target(:,:,:,ii_offset)),path_Mask_target);
            
            WriteDICOMs(B(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
            WriteDICOMs(PCi(:,:,:,ii_offset),Slices,dicomHeader,path_target_DICOM);

            [outcome_code, outcome_log] = system(['"' filepath '\MITK\matchR.exe'  '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm'  '" "' path_target_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_RegistrationAlgorhithm '" -o "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_Mask '" -m "' path_Mask_target '"']);
            if outcome_code ~= 0
                disp(outcome_log)
                warning('Error occured during image registration. See log file above for more information.\n')
            end
            T(:,:,ii_offset) = affine3d(T_conversion).invert.T*Read_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'])*T_conversion;    
        end
        
        for ii_offset = size(B,4):-1:1
            T(:,:,ii_offset) = T(:,:,ii_offset)*affine3d(T(:,:,1)).invert.T;
            WriteDICOMs(B(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
            Write_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'],T_conversion*T(:,:,ii_offset)*affine3d(T_conversion).invert.T)
            [~, ~] = system(['"'  filepath '\MITK\mapR.exe' '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_target_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" -o "' path_MoCo_DICOM 'RegImage.nrrd' '"']);
            C(:,:,Slices,ii_offset) = nrrdread([path_MoCo_DICOM 'RegImage.nrrd']);
        end
        
    end

    if d_PCA > 0
        for ii_offset = 1:size(B,4)
            T(:,:,ii_offset) = T_coarse(:,:,ii_offset)*T(:,:,ii_offset);
        end
    end
    
    % Clean up and remove temporary directories
    rmdir(path_moving_DICOM,'s')
    rmdir(path_MoCo_DICOM,'s')
    rmdir(path_target_DICOM,'s')
    if d_PCA >0
        delete(path_Mask)
        delete(path_Mask_target)
    end
    
end