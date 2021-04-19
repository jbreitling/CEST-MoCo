function T =  Read_CoRegParameter(RegFilePath)
    % Read MITK CoReg parameters from mapR (XML) file 

    % Load MAPR file as XML
    DOMnode = xmlread(RegFilePath);
    Kernels = DOMnode.getElementsByTagName('Kernel');

    %Create Transformation matrix with MITK notation: y = Tx
    T = NaN(4);
    T(4,:) = [0 0 0 1];

    for ii_direction = 1:2 %forward and inverse transformation
        Kernel = Kernels.item(ii_direction-1);

        % Inverse transformation is relevant as it defines the motion
        if strcmpi(Kernel.getAttribute('ID'), 'inverse')

            % Read rotation from the Matrix
            MatrixList = Kernel.getElementsByTagName('Matrix');
            Matrix = MatrixList.item(0);
            ValueList = Matrix.getElementsByTagName('Value');

            for ii_rotmat =1:9
                Value = ValueList.item(ii_rotmat-1);
                Column_index = str2double(Value.getAttribute('Column'))+1;
                Row_index = str2double(Value.getAttribute('Row'))+1;
                Text = Value.item(0);
                T(Row_index,Column_index) = str2double(Text.getTextContent);
            end

            % Read translation from the offset
            OffsetList = Kernel.getElementsByTagName('Offset');
            Offset = OffsetList.item(0);
            ValueList = Offset.getElementsByTagName('Value');

            for ii_transvec =1:3
                Value = ValueList.item(ii_transvec-1);
                Row_index = str2double(Value.getAttribute('Row'))+1;
                Text = Value.item(0);
                T(Row_index,4) = str2double(Text.getTextContent);
            end
        end
    end

    %Change to Matlab notation: y = xT
    T = T';
end