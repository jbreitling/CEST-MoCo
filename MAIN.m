%% Parameters for the simulation:

NumberOfRepetitions = 100;
bool_SuddenMovement = true;
bool_SubjectSpecificMovement = true;

std_rot = 0.25; % [Â°]
std_trans = 0.25; % [mm]


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(fileparts(which(mfilename)))
addpath('.\General functionality')
addpath('.\Motion corrections')

load('.\Data\GroundTruth_Info.mat')

NumberOfOffsets = size(UseableImage_GroundTruth, 4);

% The quantification results for the different metrics and approaches will
% be saved in a struct for each repetition
% (i.e. QuantificationResults.Approach.Metric(Repition))
QuantificationResults = struct();


% Loop over the repetetions
ii_repetition = 1;
while ii_repetition <= NumberOfRepetitions

    % load one of the ten groundtruth datasets
    load(['.\Data\GroundTruth_' num2str(mod(ii_repetition-1,10)+1) '.mat'])
    
    clearvars -except QuantificationResults M_GroundTruth Segment_GroundTruth R_spatial UseableImage_GroundTruth NumberOfRepetitions NumberOfOffsets bool_SuddenMovement bool_SubjectSpecificMovement std_rot std_trans ii_repetition

    %% Generate motion pattern (Equation 3)
    
    % (ii) 'Subject' dependent movement strength drawn from uniform
    % distribtuion (0,1), i.e. 1 = very strong movement, 0 = at rest
    if bool_SubjectSpecificMovement
        f_subj = rand(1);
    else
        f_subj = 1;
    end
    
    % (iii) Possiblity of sudden movement: 1% chance to amplify motion by
    % a factor of 10, otherwise 1 (normal movement)
    
    if bool_SuddenMovement
        f_sudden = 9.*(rand(NumberOfOffsets, 1)> 0.99)+ones(NumberOfOffsets,1);
    else
        f_sudden = ones(NumberOfOffsets,1);
    end
    
    % (i) Perform a Gaussian random walk for all six degrees of freedom
    % with the first entry being 0
    d_x_hat = cumsum(f_subj.*f_sudden.*std_trans.*[0; randn(NumberOfOffsets-1, 1)]);
    d_y_hat = cumsum(f_subj.*f_sudden.*std_trans.*[0; randn(NumberOfOffsets-1, 1)]);
    d_z_hat = cumsum(f_subj.*f_sudden.*std_trans.*[0; randn(NumberOfOffsets-1, 1)]);
    
    theta_x_hat = cumsum(f_subj.*f_sudden.*std_rot.*[0; randn(NumberOfOffsets-1, 1)]);
    theta_y_hat = cumsum(f_subj.*f_sudden.*std_rot.*[0; randn(NumberOfOffsets-1, 1)]);
    theta_z_hat = cumsum(f_subj.*f_sudden.*std_rot.*[0; randn(NumberOfOffsets-1, 1)]);
    
    % Calculate 4x4 transformation matrix T_hat for all offsets
    T_hat = repmat(eye(4),1,1,NumberOfOffsets); %initiliaze
    
    for ii_offset = 1:NumberOfOffsets
        %Rotation matrices R_x, R_y, R_z
        R_x = [1, 0, 0;...
            0, cosd(theta_x_hat(ii_offset)), sind(theta_x_hat(ii_offset));...
            0, -sind(theta_x_hat(ii_offset)), cosd(theta_x_hat(ii_offset))];
    
        R_y = [cosd(theta_y_hat(ii_offset)), 0, -sind(theta_y_hat(ii_offset));...
                0, 1, 0;...
                sind(theta_y_hat(ii_offset)), 0, cosd(theta_y_hat(ii_offset))];

        R_z = [cosd(theta_z_hat(ii_offset)), sind(theta_z_hat(ii_offset)), 0;...
                -sind(theta_z_hat(ii_offset)), cosd(theta_z_hat(ii_offset)), 0;...
                0, 0, 1];
            
        T_hat(1:3,1:3,ii_offset) = R_z*R_y*R_x;
        T_hat(4,1:3,ii_offset) = [d_x_hat(ii_offset), d_y_hat(ii_offset), d_z_hat(ii_offset)];
        
        
    end
    
    clear R_x R_y R_z ii_offset
    
    %% Corrupt ground truth with motion
    M_Corrupted = NaN(size(M_GroundTruth));
    Segment_Corrupted = NaN(size(Segment_GroundTruth));
    
    for ii_offset = 1:NumberOfOffsets
        M_Corrupted(:,:,:,ii_offset) = imwarp(M_GroundTruth(:,:,:,ii_offset), R_spatial, affine3d(T_hat(:,:,ii_offset)), 'FillValues', 0, 'OutputView', R_spatial);
        Segment_Corrupted(:,:,:,ii_offset) = imwarp(Segment_GroundTruth(:,:,:,ii_offset), R_spatial, affine3d(T_hat(:,:,ii_offset)), 'FillValues', 0, 'OutputView', R_spatial);
    end

    clear ii_offset
    
    % Determine 'useable', i.e. information containing part of image
    UseableImage_Corrupted = NaN(size(UseableImage_GroundTruth));
    for ii_offset = 1:NumberOfOffsets
        UseableImage_Corrupted(:,:,:,ii_offset) = imwarp(UseableImage_GroundTruth(:,:,:,ii_offset), R_spatial, affine3d(T_hat(:,:,ii_offset)), 'FillValues', 0, 'OutputView', R_spatial);
    end
    clear ii_offset
    
    UseableImage_Corrupted = ( prod(UseableImage_Corrupted,4) == 1);
    UnusableSlices = squeeze(sum(ones(size(UseableImage_Corrupted))- UseableImage_Corrupted,[1,2]) > 1000);
    UsableSlices = setdiff(1:size(M_Corrupted,3),find(UnusableSlices));
    
    % If less than 5 slices are useable, skip this repetition and start
    % over with a new motion pattern. With less than 5 slices it is not
    % reasonable to start motion correction
    if size(M_Corrupted,3) - sum(UnusableSlices) < 5
       continue
    end
    
    M_Corrupted(:,:,UnusableSlices,:) = 0;
    Segment_Corrupted(:,:,UnusableSlices,:) = 0;
    
    %% Corrupt data with Rician noise (1% of maximum image intensity)
    M_Corrupted = ricernd(M_Corrupted, 0.01*max(M_Corrupted,[],'all'));
    
    %% 
    [dicomHeader, ~] = GetDicomHeader('.\Data\Dicoms\');
    WriteDICOMs(M_Corrupted(:,:,:,1),UsableSlices,dicomHeader,'.\target_DICOM\');
    
    %% Generate reference dataset by using 'perfect' inverse transformation 
    
    %M_ref = NaN(size(M_Corrupted));
    Segment_ref = NaN(size(Segment_Corrupted));

    
    M_Corrupted = round(M_Corrupted);    
    for ii_offset = 1:NumberOfOffsets
        Segment_ref(:,:,:,ii_offset) = imwarp(Segment_Corrupted(:,:,:,ii_offset), R_spatial, invert(affine3d(T_hat(:,:,ii_offset))), 'FillValues', 0, 'OutputView', R_spatial);
    end
    
    % Generate ground truth by using 'perfect' motion correction
    M_ref = MotionCorrection_PERFECT(M_Corrupted,UsableSlices,['.\target_DICOM\Slice_' num2str(find(~UnusableSlices,1)) '.dcm'],'.\Data\Dicoms\',T_hat);
    clear ii_offset
    
    M_ref(:,:,UnusableSlices,:) = 0;
    Segment_ref(:,:,UnusableSlices,:) = 0;
    
    %% Apply different motion correction approaches
    
    % The results of the different motion correction approaches are saved
    % in a struct for easier accessibility
    MotionCorrectionResults = struct();
    
    
    try
    
    % No motion correction
    MotionCorrectionResults.Corrupted.M = M_Corrupted;
    MotionCorrectionResults.Corrupted.T = repmat(eye(4),1,1,NumberOfOffsets);
    MotionCorrectionResults.Corrupted.Time = 0;tic;
    
    % Conventional motion correction
    [MotionCorrectionResults.CONV.M,...
        MotionCorrectionResults.CONV.T] = MotionCorrection_CONV(M_Corrupted,UsableSlices,['.\target_DICOM\Slice_' num2str(find(~UnusableSlices,1)) '.dcm'],'.\Data\Dicoms\');
    MotionCorrectionResults.CONV.Time = toc;tic;
    
    % Proposed motion correction with identificication and mitigation
    [MotionCorrectionResults.PROPOSED.M,...
        MotionCorrectionResults.PROPOSED.T] = MotionCorrection_PROPOSED(M_Corrupted,UsableSlices,['.\target_DICOM\Slice_' num2str(find(~UnusableSlices,1)) '.dcm'],'.\Data\Dicoms\');
    MotionCorrectionResults.PROPOSED.Time = toc;tic;
    
    % Motion correction with linear interpolation
    [MotionCorrectionResults.LinInterp.M,...
        MotionCorrectionResults.LinInterp.T] = MotionCorrection_LinInterp(M_Corrupted,UsableSlices,['.\target_DICOM\Slice_' num2str(find(~UnusableSlices,1)) '.dcm'],'.\Data\Dicoms\',25:35);
    MotionCorrectionResults.LinInterp.Time = toc;tic;

    % LRAZ motion correction
    [MotionCorrectionResults.LRAZ.M,...
        MotionCorrectionResults.LRAZ.T,...
        MotionCorrectionResults.LRAZ.Iteration] = MotionCorrection_LRAZ(M_Corrupted,UsableSlices,'.\Data\Dicoms\',0.02,30);
    MotionCorrectionResults.LRAZ.Time = toc;tic;
    
    % RPCA+PCA_R motion correction
    [MotionCorrectionResults.RPCA_PCA_R.M,...
        MotionCorrectionResults.RPCA_PCA_R.T] = MotionCorrection_RPCA_PCA_R(M_Corrupted,UsableSlices,'.\Data\Dicoms\',0.75,0);
    MotionCorrectionResults.RPCA_PCA_R.Time = toc;
    

    catch
        disp('Error')
        continue
    end
    
    %% Visualize results
    figure, hold on
    subplot(3,3,1),hold on, plot(d_x_hat,'Color','k'), box on, ylabel('{\Delta}x [mm]')
    subplot(3,3,2),hold on, plot(d_y_hat,'Color','k'), box on, ylabel('{\Delta}y [mm]')
    subplot(3,3,3),hold on, plot(d_z_hat,'Color','k'), box on, ylabel('{\Delta}z [mm]')
    subplot(3,3,4),hold on, plot(theta_x_hat,'Color','k'), box on, ylabel('Pitch [°]')
    subplot(3,3,5),hold on, plot(theta_y_hat,'Color','k'), box on, ylabel('Roll [°]')
    subplot(3,3,6),hold on, plot(theta_z_hat,'Color','k'), box on, ylabel('Yaw [°]')
    subplot(3,3,6),hold on, plot(theta_z_hat,'Color','k'), box on, ylabel('Yaw [°]')
    subplot(3,3,9),hold on, plot(1,1,'Color','k'), box on, xticks([]), yticks([])

    Approaches = fieldnames(MotionCorrectionResults);
    colors = lines(numel(Approaches));
    for ii_approach = 1:numel(Approaches)
        d_x = NaN(NumberOfOffsets,1);d_y = d_x; d_z = d_x; theta_x = d_x; theta_y = d_x; theta_z = d_x;
        
        for ii_offset = 1:NumberOfOffsets
            [d_x(ii_offset),d_y(ii_offset),d_z(ii_offset),...
                theta_x(ii_offset),theta_y(ii_offset),theta_z(ii_offset)] =...
                calculateTransformationParameters(MotionCorrectionResults.(Approaches{ii_approach}).T(:,:,ii_offset));
        end
        
        subplot(3,3,1), plot(d_x,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
        subplot(3,3,2), plot(d_y,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
        subplot(3,3,3), plot(d_z,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
        subplot(3,3,4), plot(theta_x,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
        subplot(3,3,5), plot(theta_y,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
        subplot(3,3,6), plot(theta_z,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
        subplot(3,3,9), plot(1,1,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets])
    end
    legend([{'Applied'};Approaches])
    
    clear ii_approach ii_offset d_x d_y d_z theta_x theta_y theta_z
    
    %% Transformation of the segments
    
    Approaches = fieldnames(MotionCorrectionResults);
    
    for ii_approach = 1:numel(Approaches)
        Segment = NaN(size(Segment_Corrupted));
        
        for ii_offset = 1:NumberOfOffsets
            Segment(:,:,:,ii_offset) = imwarp(Segment_Corrupted(:,:,:,ii_offset), R_spatial, invert(affine3d(MotionCorrectionResults.(Approaches{ii_approach}).T(:,:,ii_offset))), 'FillValues', 0, 'OutputView', R_spatial);
        end
        
        Segment = (Segment == 1);
        Segment(:,:,UnusableSlices,:) = 0;
        
        MotionCorrectionResults.(Approaches{ii_approach}).Segment = Segment;
    end
    
    clear Segment ii_approach ii_offset

    
    %% Quantification of the motion correction
    
    for ii_approach = 1:numel(Approaches)
        
        % Spectral error - normalized root-mean-squared error (NRMSE)
        NRMSE_i = NaN(1,NumberOfOffsets);
        
        for ii_offset = 1:NumberOfOffsets
            Segment_intersection = (MotionCorrectionResults.(Approaches{ii_approach}).Segment + Segment_ref == 2);
            Segment_intersection(:,:,:,setdiff(1:NumberOfOffsets,ii_offset)) = false;
            NRMSE_i(ii_offset) = sqrt(mean(((MotionCorrectionResults.(Approaches{ii_approach}).M(Segment_intersection) - M_ref(Segment_intersection))./M_ref(Segment_intersection(:,:,:,[ii_offset setdiff(1:NumberOfOffsets,ii_offset)]))).^2));
        end

        NRMSE = sqrt(mean(NRMSE_i.^2));
        
        % Image misalignment - maximum of root-mean-squared deviation (d_RMS)
        R = 70; %[mm]
        d_RMS_i = NaN(1,NumberOfOffsets);
        
        for ii_offset = 1:NumberOfOffsets
            M = MotionCorrectionResults.(Approaches{ii_approach}).T(:,:,ii_offset)/T_hat(:,:,ii_offset) - eye(4);
            d_RMS_i(ii_offset) = sqrt(1/5 * R^2 * trace((M(1:3,1:3))'*M(1:3,1:3)) + (M(4,1:3))*(M(4,1:3))');
        end

        d_RMS = max(d_RMS_i);
        
        % Append the quantification results to the previous results
        if ~isfield(QuantificationResults,Approaches{ii_approach})
            QuantificationResults.(Approaches{ii_approach}).NRMSE = NRMSE;
            QuantificationResults.(Approaches{ii_approach}).d_RMS = d_RMS;
            QuantificationResults.(Approaches{ii_approach}).Time = MotionCorrectionResults.(Approaches{ii_approach}).Time;
        else
            QuantificationResults.(Approaches{ii_approach}).NRMSE(end+1) = NRMSE;
            QuantificationResults.(Approaches{ii_approach}).d_RMS(end+1) = d_RMS;
            QuantificationResults.(Approaches{ii_approach}).Time(end+1) = MotionCorrectionResults.(Approaches{ii_approach}).Time;
        end
        
        % Plot in addition results of quantification
        subplot(3,3,7),hold on, plot(d_RMS_i,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets]), box on, ylabel('d_{RMS} [mm]')
        subplot(3,3,8),hold on, plot(NRMSE_i*100,'Color',colors(ii_approach,:));xlim([1 NumberOfOffsets]), box on, ylabel('NRMSE_{i} [%]')
        
        if isfield(MotionCorrectionResults.(Approaches{ii_approach}), 'Iteration')
            if ~isfield(QuantificationResults.(Approaches{ii_approach}),'Iteration')
                QuantificationResults.(Approaches{ii_approach}).Iteration = MotionCorrectionResults.(Approaches{ii_approach}).Iteration;
            else
                QuantificationResults.(Approaches{ii_approach}).Iteration(end+1) = MotionCorrectionResults.(Approaches{ii_approach}).Iteration;
            end
        end
        
    end
    hold off
    clear ii_approach ii_offset Segment_intersection_i M_CoReg_i M_ref_i NRMSE d_RMS R M NRMSE_i d_RMS_i colors

    
    %% Increase iterator
    ii_repetition = ii_repetition +1;
end



%% Plot the results of the simulation

Approaches = fieldnames(QuantificationResults);
Metrics = fieldnames(QuantificationResults.(Approaches{1}));

colors = lines(numel(Approaches));

for ii_metric = 1:numel(Metrics)
    figure, hold on

    plotData = NaN(NumberOfRepetitions,numel(Approaches));        
    for ii_approach = 1:numel(Approaches)
        temp = QuantificationResults.(Approaches{ii_approach}).(Metrics{ii_metric});
        plotData(1:(numel(temp)),ii_approach) = temp;
    end

    boxplot(plotData,'Whisker',100,'notch','off','Colors',colors,'Widths',0.7,'jitter',0.3,'Symbol','k+');
    
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(end+1-j,:),'FaceAlpha',.5);
        set(h(j,:),'LineWidth',1);
        set(h(j,:),'LineStyle','-');
    end
    clear h j

    hold off
end

clear Approaches Metrics colors ii_metric plotData ii_approach
