function nrrdwrite(filename, X, meta, nrrdVersion)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NRRDWRITE  Export NRRD imagery and metadata.
% NRRDWRITE(FILENAME, X, META, nrrdVersion) writes the image volume X and associated
% metadata to a file specified by FILENAME in the NRRD-format.
%
% Compressed NRRD is not supported at the moment (will be written as RAW)!
%
% See the format specification online:
% http://teem.sourceforge.net/nrrd/format.html
%
% Copyright 2012 Dr. Sarah Mang
%
%
% SEEALSO nrrdread
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open file.
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

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

% Check if required fields are provided
assert(isfield(meta, 'sizes') && ...
       isfield(meta, 'dimension') && ...
       isfield(meta, 'encoding') && ...
       isfield(meta, 'endian'), ...
       'Missing required metadata fields.')

% get back from MATLAB's order.
%if(size(X) ~= meta.sizes) 
 %   X = permute(X, [2 1 3 4]);
%end
% print header
% write as text
fprintf(fid,[nrrdVersion,'\n']);
if length(size(X)) >=3
    meta.sizes = num2str(size(X,1));
    for i = 2: length(size(X))
        meta.sizes = [meta.sizes,' ',num2str(size(X,i))];
    end
    meta.dimension = num2str(length(size(X)));
else
    meta.sizes = num2str(size(X,1));
    for i = 2: length(size(X))
        meta.sizes = [meta.sizes,' ',num2str(size(X,i))];
    end
    for i = length(size(X))+1:3
        meta.sizes = [meta.sizes,' 1'];
    end
    meta.dimension = '3';
end

% ziped bodies are not supported at the moment
if strcmp(meta.encoding, 'gzip') | strcmp(meta.encoding, 'gz')
    meta.encoding = 'raw';
end
% write meta data
fn = fieldnames(meta);
for i = 1:max(size(fn))
    value = getfield(meta,fn{i});
    assert(isstr(value), 'Meta Data needs to be all strings')
    if(strcmp(fn{i}, 'org_mitk_timegeometry_timepoints') || strcmp(fn{i}, 'org_mitk_timegeometry_type'))
        metaLine = [strtrim(strrep(fn{i}, 'space', 'space ')),':=',value,'\n'];
    else
        metaLine = [strtrim(strrep(fn{i}, 'space', 'space ')),': ',value,'\n'];
    end
    fprintf(fid,metaLine);
end
fprintf(fid,'\n');   

% print data
% write in binary format
data = reshape(X, 1, []); %reshape for better performance
data = adjustEndian(data, meta);
writeData(fid, data, meta);


function writeData(fidIn, data, meta)
datatype=meta.type;
switch (meta.encoding)
 case {'raw'}
  
  ok = fwrite(fidIn, data, meta.type);
  
 case {'gzip'}
     
     % Store in a raw file before compressing
     tmpBase = tempname(pwd);
     tmpFile = [tmpBase '.gz'];
     fidTmpRaw = fopen(tmpBase, 'wb');
     assert(fidTmpRaw > 3, 'Could not open temporary file for GZIP compression');
     
     fwrite(fidTmpRaw, data(:), datatype);
     fclose(fidTmpRaw);
     
     % Now we gzip our raw file
     gzip(tmpBase);
     
     % Finally, we put this info into our nrrd file (fidIn)
     fidTmpRaw = fopen(tmpFile, 'rb');
     tmp = fread(fidTmpRaw, inf, [datatype '=>' datatype]);
     cleaner = onCleanup(@() fclose(fidTmpRaw));
     ok = fwrite (fidIn, tmp, datatype);
     
     delete (tmpBase);
     delete (tmpFile);


 case {'txt', 'text', 'ascii'}
  
  ok = fwrite(fidIn, data, meta.type);
  
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

