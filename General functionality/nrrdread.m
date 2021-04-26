function [X, meta, nrrdVersion] = nrrdread(filename)
%NRRDREAD  Import NRRD imagery and metadata.
%   [X, META] = NRRDREAD(FILENAME) reads the image volume and associated
%   metadata from the NRRD-format file specified by FILENAME.
%
%   Current limitations/caveats:
%   * "Block" datatype is not supported.
%   * Only tested with "gzip" and "raw" file encodings.
%   * Very limited testing on actual files.
%   * I only spent a couple minutes reading the NRRD spec.
%
%   See the format specification online:
%   http://teem.sourceforge.net/nrrd/format.html
%
% Copyright 2012 The MathWorks, Inc.


% Open file.
fid = fopen(filename, 'rb');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

% Magic line.
theLine = fgetl(fid);
assert(numel(theLine) >= 4, 'Bad signature in file.')
assert(isequal(theLine(1:4), 'NRRD'), 'Bad signature in file.')
nrrdVersion = theLine;

% The general format of a NRRD file (with attached header) is:
% 
%     NRRD000X
%     <field>: <desc>
%     <field>: <desc>
%     # <comment>
%     ...
%     <field>: <desc>
%     <key>:=<value>
%     <key>:=<value>
%     <key>:=<value>
%     # <comment>
% 
%     <data><data><data><data><data><data>...

meta = struct([]);

% Parse the file a line at a time.
while (true)

  theLine = fgetl(fid);
  
  if (isempty(theLine))
    % End of the header.
    break;
  end
  
  if (isequal(theLine(1), '#'))
      % Comment line.
      continue;
  end
  
  % "fieldname:= value" or "fieldname: value" or "fieldname:value"
  parsedLine = regexp(theLine, ':=?\s*', 'split','once');
  
  assert(numel(parsedLine) == 2, 'Parsing error')
  
  field = lower(parsedLine{1});
  value = parsedLine{2};
  
  field(isspace(field)) = '';
  field(find(field=='-'))=[]; %eliminate "-", which is not allowed in field names
  meta(1).(field) = value;
  
  if (feof(fid))
    % End of the header.
    break;
  end
  
end

datatype = getDatatype(meta.type);

% Get the size of the data.
% assert(isfield(meta, 'sizes') && ...
%        isfield(meta, 'dimension') && ...
%        isfield(meta, 'encoding') && ...
%        isfield(meta, 'endian'), ...
%        'Missing required metadata fields.')
% MITK does not always include endianess 
assert(isfield(meta, 'sizes') && ...
       isfield(meta, 'dimension') && ...
       isfield(meta, 'encoding') , ...
       'Missing required metadata fields.')
if ~isfield(meta,'endian')
    % todo need to check what MITK assumes
    meta.endian = 'big';
end

dims = sscanf(meta.sizes, '%d');
ndims = sscanf(meta.dimension, '%d');
assert(numel(dims) == ndims);

[pathstr, name, ext] = fileparts(filename); 
data = readData(fid, meta, datatype, pathstr);
data = adjustEndian(data, meta);

% Reshape to matrix.
X = reshape(data, dims');

% If DWI: Check Ordering
if isfield(meta, 'modality') && strcmp( meta.modality,'DWMRI')
    dimType = textscan(meta.kinds, '%s');
    spaceDir =  textscan(meta.spacedirections, '%s');
    % find gradient domain:
    gradDim = 0;
    for i =1:numel(dims)
        if strcmp(dimType{1}{i},'list') | strcmp(dimType{1}{i},'vector')
            gradDim = i;
            break
        end
    end
    assert(gradDim > 0, 'Gradient dimension could not be identified from header')
    assert(numel(dims)<5 && numel(dims)>2, 'unsuported number of dimensions for DWI (needs to be 3 or 4)')
    % desired ordering is [x,y,z,g]
    if numel(dims) == 4
        switch gradDim
            case 1
                permOrder = [3 2 4 1];
                spaceOrder= [2 3 4 1];
            case 2
                permOrder = [3 1 4 2];
                spaceOrder= [1 3 4 2];
            case 3
                permOrder = [2 1 4 3];
                spaceOrder= [1 2 4 3];
            otherwise
                permOrder = [2 1 3 4];
                spaceOrder= [1 2 3 4];
        end
        dimType = dimType{1}(permOrder);
        spaceDir = spaceDir{1}(spaceOrder);
        meta.kinds = [dimType{1},' ',dimType{2},' ',dimType{3},' ',dimType{4}];
        meta.spacedirections = [spaceDir{1},' ',spaceDir{2},' ',spaceDir{3},' ',spaceDir{4}];
    else
        switch gradDim
            case 1
                permOrder = [3 2 1]; 
                spaceOrder= [2 3 1];
            case 2
                permOrder = [3 1 2];
                spaceOrder= [1 3 2];
            otherwise
                permOrder = [2 1 3];
                spaceOrder= [1 2 3];
        end
        dimType = dimType{1}(permOrder);
        spaceDir = spaceDir{1}([1,2,3]);
        spaceDir = spaceDir{1}(spaceOrder);
        meta.spacedirections = [spaceDir{1},' ',spaceDir{2},' ',spaceDir{3},' ',spaceDir{4}];
    end
    X = permute(X, permOrder);
else
    % Get into MATLAB's order.
    X = permute(X, [2 1 3]);
end


function datatype = getDatatype(metaType)

% Determine the datatype
switch (metaType)
 case {'signed char', 'int8', 'int8_t'}
  datatype = 'int8';
  
 case {'uchar', 'unsigned char', 'uint8', 'uint8_t'}
  datatype = 'uint8';

 case {'short', 'short int', 'signed short', 'signed short int', ...
       'int16', 'int16_t'}
  datatype = 'int16';
  
 case {'ushort', 'unsigned short', 'unsigned short int', 'uint16', ...
       'uint16_t'}
  datatype = 'uint16';
  
 case {'int', 'signed int', 'int32', 'int32_t'}
  datatype = 'int32';
  
 case {'uint', 'unsigned int', 'uint32', 'uint32_t'}
  datatype = 'uint32';
  
 case {'longlong', 'long long', 'long long int', 'signed long long', ...
       'signed long long int', 'int64', 'int64_t'}
  datatype = 'int64';
  
 case {'ulonglong', 'unsigned long long', 'unsigned long long int', ...
       'uint64', 'uint64_t'}
  datatype = 'uint64';
  
 case {'float'}
  datatype = 'single';
  
 case {'double'}
  datatype = 'double';
  
 otherwise
  assert(false, 'Unknown datatype')
end



function data = readData(fidIn, meta, datatype, pathstr)
% required if header is in separate file
if isfield(meta, 'datafile') && length(meta.datafile) > 0
    fidIn = fopen(fullfile(pathstr,meta.datafile));
    cleaner2 = onCleanup(@() fclose(fidIn));
end

switch (meta.encoding)
 case {'raw'}
  
  data = fread(fidIn, inf, [datatype '=>' datatype]);
  
 case {'gzip', 'gz'}

  % required if header is in separate file
  tmpBase = tempname();
  tmpFile = [tmpBase '.gz'];
  if isfield(meta, 'datafile') && length(meta.datafile) > 0
      copyfile(fullfile(pathstr,meta.datafile),tmpFile)
  
      meta.datafile = '';
  else
      fidTmp = fopen(tmpFile, 'wb');
      assert(fidTmp > 3, 'Could not open temporary file for GZIP decompression')
      
      tmp = fread(fidIn, inf, 'uint8=>uint8');
      fwrite(fidTmp, tmp, 'uint8');
      fclose(fidTmp);
  end
  
  gunzip(tmpFile)
      
  fidTmp = fopen(tmpBase, 'rb');
  cleaner = onCleanup(@() fclose(fidTmp));
  
  meta.encoding = 'raw';
  data = readData(fidTmp, meta, datatype, pathstr);
  
  
 case {'txt', 'text', 'ascii'}
  
  data = fscanf(fidIn, '%f');
  data = cast(data, datatype);
  
 otherwise
  assert(false, 'Unsupported encoding')
end



function data = adjustEndian(data, meta)

[~,~,endian] = computer();

needToSwap = (isequal(endian, 'B') && isequal(lower(meta.endian), 'little')) || ...
             (isequal(endian, 'L') && isequal(lower(meta.endian), 'big'));
         
if (needToSwap)
    data = swapbytes(data);
end