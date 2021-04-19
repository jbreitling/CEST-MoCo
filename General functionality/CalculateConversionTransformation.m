function T_conversion = CalculateConversionTransformation(dicomHeader)

    % Calculate T_rot and T_shift, which yield transformation from
    % DICOM reference system to the image system, i.e. centered at the
    % middle of the image

    T_rot = diag(ones(4,1));
    T_rot(1:3,1) = dicomHeader{1, 1}.ImageOrientationPatient(1:3);
    T_rot(1:3,2) = dicomHeader{1, 1}.ImageOrientationPatient(4:6);
    T_rot(1:3,3) = cross(T_rot(1:3,1),T_rot(1:3,2));

    T_shift = diag(ones(4,1));

    shift_to_center = zeros(3,1);
    shift_to_center(2) = (double(dicomHeader{1, 1}.Rows)/2-0.5) *dicomHeader{1, 1}.PixelSpacing(1);
    shift_to_center(1) = (double(dicomHeader{1, 1}.Columns)/2 -0.5) *dicomHeader{1, 1}.PixelSpacing(1);
    shift_to_center = T_rot(1:3,1:3)*shift_to_center;

    z = NaN(numel(dicomHeader),1); y = z; x = z;

    for ii_files = 1: numel(dicomHeader)
        z(ii_files) = dicomHeader{1, ii_files}.ImagePositionPatient(3);
        y(ii_files) = dicomHeader{1, ii_files}.ImagePositionPatient(2);
        x(ii_files) = dicomHeader{1, ii_files}.ImagePositionPatient(1);
    end

    center = [ mean(x); mean(y); mean(z)];
    T_shift(4,1:3) = - center - shift_to_center;
    
    T_conversion = T_shift*T_rot;
end