function Write_CoRegParameter(RegFilePath,T)
% Function which allows to write new CoReg Parameters to an existing mapR
% file

    % Change to MITK notation: y = Tx
    T_inv = T';
    InverseRotMat = T_inv(1:3,1:3);
    InverseTransVec = T_inv(1:3,4);

    T_forw = affine3d(T).invert.T';
    ForwardRotMat = T_forw(1:3,1:3);
    ForwardTransVec = T_forw(1:3,4);

    % Load MAPR file as XML
    DOMnode = xmlread(RegFilePath);
    Kernels = DOMnode.getElementsByTagName('Kernel');

    for ii_direction = 1:2 %forward and inverse transformation
        Kernel = Kernels.item(ii_direction-1);

        if strcmpi(Kernel.getAttribute('ID'), 'direct') %forward transformation
            RotMat = ForwardRotMat;
            TransVec = ForwardTransVec;
        elseif strcmpi(Kernel.getAttribute('ID'), 'inverse') %inverse transformation
            RotMat = InverseRotMat;
            TransVec = InverseTransVec;
        end

        % Save rotation to Matrix
        MatrixList = Kernel.getElementsByTagName('Matrix');
        Matrix = MatrixList.item(0);
        ValueList = Matrix.getElementsByTagName('Value');

        for ii_rotmat =1:9
            Value = ValueList.item(ii_rotmat-1);
            Column_index = str2double(Value.getAttribute('Column'))+1;
            Row_index = str2double(Value.getAttribute('Row'))+1;
            Text = Value.item(0);
            Text.setTextContent(num2str(RotMat(Row_index,Column_index),'%.12f'));
        end

        % Save rotation to Matrix string
        MatrixStrList = Kernel.getElementsByTagName('MatrixStr');
        MatrixStr = MatrixStrList.item(0);
        Text = MatrixStr.item(0);
        Text.setTextContent(strjoin(string((RotMat'))))

        % Save translation to offset
        OffsetList = Kernel.getElementsByTagName('Offset');
        Offset = OffsetList.item(0);
        ValueList = Offset.getElementsByTagName('Value');

        for ii_transvec =1:3
            Value = ValueList.item(ii_transvec-1);
            Row_index = str2double(Value.getAttribute('Row'))+1;
            Text = Value.item(0);
            Text.setTextContent(num2str(TransVec(Row_index),'%.12f'));
        end

        % Save translation as string
        OffsetStrList = Kernel.getElementsByTagName('OffsetStr');
        OffsetStr = OffsetStrList.item(0);
        Text = OffsetStr.item(0);
        Text.setTextContent(strjoin(string((TransVec(:)))))
    end

    % Replace the RegFile
    xmlwrite(RegFilePath,DOMnode)
end