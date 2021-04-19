function [M_CoReg, T, ii_iteration] = MotionCorrection_LRAZ(M_Corrupted,Slices,path_original_DICOM,tau,NumberOfIterations)
%MotionCorrection_LRAZ  Cutting-edge motion correction exploiting a
% low-rank approximation of the z-spectrum (LRAZ)
%   Cutting-edge motion correction exploiting a low-rank approximation of
%   the z-spectrum (LRAZ) as proposed in:
%
%   'Wech T, Köstler H. Robust Motion Correction In CEST Imaging Exploiting
%   Low-Rank Approximation Of The Z-Spectrum.
%   Magn Reson Med. 2018;80(5):1979-1988.'
%
%   The variable nominclature follows the original paper.The algorithm 
%   co-registers the motion corrupted images M_Corrupted (x,y,z,dw)to each
%   other using only the slices specified . As MITK requires DICOMs the
%   folder of an original DICOM is required: path_original_DICOM.
%   Returns the co-registered data M_CoReg and the 4x4 transformation
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
    
    if exist(path_MoCo_DICOM,'dir')
        rmdir(path_MoCo_DICOM,'s')
    end
    mkdir(path_MoCo_DICOM);
    
    T = repmat(eye(4),1,1,size(M_Corrupted,4));

    R_spatial = imref3d(size(M_Corrupted(:,:,:,1)),...
        size(M_Corrupted,2)*dicomHeader{1}.PixelSpacing(2).*[-0.5 +0.5],...
        size(M_Corrupted,1)*dicomHeader{1}.PixelSpacing(1).*[-0.5 +0.5],...
        size(M_Corrupted,3)*dicomHeader{1}.SliceThickness.*[-0.5 +0.5]);

    Mask_orig = zeros(size(M_Corrupted(:,:,:,1)));
    Mask_orig(:,:,Slices) = 1;
    
    tau_ii = interp1([1 NumberOfIterations],[tau 0],1:NumberOfIterations, 'linear');
    Mask_offset = repmat(Mask_orig,1,1,1,size(M_Corrupted,4));
    
    M_CoReg = M_Corrupted;
    
    for ii_iteration = 1:NumberOfIterations
        
        for ii_offset =  1:size(M_Corrupted,4)
            Mask_offset(:,:,:,ii_offset) = imwarp(Mask_orig, R_spatial, invert(affine3d(T(:,:,ii_offset))), 'FillValues', 0, 'OutputView', R_spatial);
        end
        
        Mask_ii = (prod(Mask_offset,4) == 1);
        
        if sum(Mask_ii,'all') == 0
            warning(['Mask of LRAZ algorithm is empty. No more overlap between different offsets. Algorithm was stopped after ' num2str(ii_iteration-1) ' iterations!'])
            break
        end
        
        WriteMask(dicomHeader,double(Mask_ii),path_Mask);
        
        M_ii = permute(M_CoReg,[4 1 2 3]);
        C_ii = M_ii(:,Mask_ii)';
        [U_ii,S_ii,V_ii] = svd(C_ii,'econ');
        
        S_ii(S_ii < tau_ii(ii_iteration) .* S_ii(1,1)) = 0;
        S_ii(S_ii ~=0) = S_ii(S_ii ~=0) - tau_ii(ii_iteration)* S_ii(1,1);
        
        L_ii = U_ii*S_ii*V_ii';
        
        L_ii_z = zeros(size(M_ii));
        L_ii_z(:,Mask_ii) = L_ii';
        L_ii_z = permute(L_ii_z,[2 3 4 1]);
        
        for ii_offset = 1:size(M_Corrupted,4)
            WriteDICOMs(M_Corrupted(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
            WriteDICOMs(L_ii_z(:,:,:,ii_offset),Slices,dicomHeader,path_target_DICOM);

            % Matching
            [outcome_code, outcome_log] = system(['"' filepath '\MITK\matchR.exe'  '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm'  '" "' path_target_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_RegistrationAlgorhithm '" -o "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_Mask '"']);
            if outcome_code ~= 0
                disp(outcome_log)
                warning('Error occured during image registration. See log file above for more information.\n')
                % Clean up and remove temporary directories
                rmdir(path_moving_DICOM,'s')
                rmdir(path_MoCo_DICOM,'s')
                rmdir(path_target_DICOM,'s')
                delete(path_Mask)
                break
            end
            
            % Read transformation matrix
            T(:,:,ii_offset) = affine3d(T_conversion).invert.T*Read_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'])*T_conversion;    
        end
        
        if outcome_code ~= 0
            break
        end
        
        % Mapping
        for ii_offset = size(M_Corrupted,4):-1:1
            T(:,:,ii_offset) = T(:,:,ii_offset)*affine3d(T(:,:,1)).invert.T;
            WriteDICOMs(M_Corrupted(:,:,:,ii_offset),Slices,dicomHeader,path_moving_DICOM);
            Write_CoRegParameter([path_MoCo_DICOM 'RegFile.mapr'],T_conversion*T(:,:,ii_offset)*affine3d(T_conversion).invert.T)
            [~, ~] = system(['"'  filepath '\MITK\mapR.exe' '" "' path_moving_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" "' path_MoCo_DICOM 'RegFile.mapr' '" -t "' path_target_DICOM 'Slice_' num2str(Slices(1)) '.dcm' '" -o "' path_MoCo_DICOM 'RegImage.nrrd' '"']);
            M_CoReg(:,:,Slices,ii_offset) = nrrdread([path_MoCo_DICOM 'RegImage.nrrd']);
        end
    end

    % Clean up and remove temporary directories
    rmdir(path_moving_DICOM,'s')
    rmdir(path_MoCo_DICOM,'s')
    rmdir(path_target_DICOM,'s')
    delete(path_Mask)
    
end