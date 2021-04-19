function WriteMask(dicomHeader,Mask,Path_Mask)
% Function to write Segment/Mask as a NRRD-file. Can be used to create the
% input for the MITK image registration.
    
    if iscell(dicomHeader)
        for ii_cell = 1:numel(dicomHeader)
            if dicomHeader{ii_cell}.AcquisitionNumber == 1 && dicomHeader{ii_cell}.InstanceNumber == 1
                dicomHeader = dicomHeader{ii_cell};
                break
            end
        end
    end

    % orientation of image - nrrd scales the vectors by the pixel
    % spacing
    directions(1:3,1) = dicomHeader.ImageOrientationPatient(1:3);
    directions(1:3,2) = dicomHeader.ImageOrientationPatient(4:6);
    directions(1:3,3) = cross(directions(1:3,1),directions(1:3,2));
    directions = directions.*[dicomHeader.PixelSpacing(1), dicomHeader.PixelSpacing(2), dicomHeader.SliceThickness];

    % define the origin of the image
    origin = dicomHeader.ImagePositionPatient;

    %create nrrd header
    meta.dimension = '3';
    meta.space = 'left-posterior-superior';
    meta.endian = 'little';
    meta.encoding = 'gzip';
    meta.kinds = 'domain domain domain';
    meta.type ='unsigned short';
    meta.sizes = '0 0 0';
    meta.spacedirections = append('(',num2str(directions(1,1)),',',num2str(directions(2,1)),',',num2str(directions(3,1)),') (',num2str(directions(1,2)),',',num2str(directions(2,2)),',',num2str(directions(3,2)),') (',num2str(directions(1,3)),',',num2str(directions(2,3)),',',num2str(directions(3,3)),')');
    meta.spaceorigin = append('(',num2str(origin(1)),',',num2str(origin(2)),',',num2str(origin(3)),')');

    % Image has to be reoriented
    Mask = flip(rot90(Mask,3),2);

    % write nrrd file
    nrrdVersion = 'NRRD0004';
    nrrdwrite(Path_Mask,Mask,meta,nrrdVersion);   
end