function [d_x,d_y,d_z,theta_x,theta_y,theta_z] = calculateTransformationParameters(T)
%Author: JB (based on PBD work)
%Function to calculate the rotation angles and translations from the 4x4 matrix
%IMPORTANT: d_x, d_y and d_z are defined in millimeters! Thetas in deg.
%IMPORTANT: T has to be defined as a right-handed matrix, i.e. in the
% MatLab style y = xT 


        % check that matrix is of size 4x4
        if (~isequal(size(T), [4 4]))
            error('Error');
        end
        
        %Copy of PBD calculations from imageRegistration
        RotMat(1:3,1:3) = T(1:3,1:3)';
        TransVec =T(4,1:3);

        % calculate euler angles
        %euler_radian = rotm2eul(RotMat,'ZYX'); % euler angles in radian 'ZYX'
        sy = sqrt(RotMat(1,1)^2 + RotMat(2,1)^2);
        if(~(sy < 1e-6))
            theta_x = atan2d(RotMat(3,2),RotMat(3,3));  % rotation around x
            theta_y = atan2d(-RotMat(3,1),sy);          % rotation around y
            theta_z = atan2d(RotMat(2,1),RotMat(1,1));  % rotation around z
        else
            theta_x = atan2d(-RotMat(2,3),RotMat(2,2)); % rotation around x
            theta_y = atan2d(-RotMat(3,1),sy);          % rotation around y
            theta_z = 0;                                % rotation around z
        end

        % Translations:
        d_x = TransVec(1); % translation in x [mm]
        d_y = TransVec(2); % translation in y [mm]
        d_z = TransVec(3); % translation in z [mm]

end