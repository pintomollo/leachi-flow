function [nframes, ssize] = size_data(fname)
% SIZE_DATA extracts the number of frames as well as the size of a frame from a file.
%
%   [NFRAMES, IMG_SIZE] = SIZE_DATA(FNAME) returns the NFRAMES and the IMG_SIZE from
%   FNAME. It requires a file format readable by imformats. SIZE_DATA returns -1 in
%   case of error.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  % Initialize the outputs
  nframes = -1;
  ssize = NaN(1,2);

  % If no filename is provided, we have a problem
  if (isempty(fname))
    warning('Tracking:SizeData', 'No file name was provided.');
    return;

  % Or, potentially is it's a structure
  elseif(isstruct(fname))

    % We're looking only for one field with the name 'fname', everyhting else if wrong
    if (~isfield(fname, 'fname') || isempty(fname.fname))
      warning('Tracking:SizeData', 'No file name was found in the provided structure.');
      return;

    % Try to open it
    else
      fid = fopen(fname.fname, 'r');
    end

  % If we have a string, try to open the file
  elseif (ischar(fname))
    fid = fopen(fname, 'r');

  % Otherwise, we do not know what to do
  else
    warning('Tracking:SizeData', 'Unable to extract a filename from an "%s" object.', class(fname));
    return
  end

  % Here, we could not open the file previously, so problem...
  if (fid == -1)
    warning('Tracking:SizeData', 'Unable to open file "%s" for reading.', fname);
    return

  % Now we can work ! Most of this part was extracted from the code of imfinfo
  % for a maximal speedup (I got rid of quite some checking)
  else
    % Get the full filename just in case
    filename = fopen(fid);
    fclose(fid);

    % Get ready to extract the image format
    format = '';

    % Look for the extension
    idx = find(filename == '.');
    if (~isempty(idx))

      % Extract the extension
      extension = lower(filename(idx(end)+1:end));

      % Look up the extension in the file format registry.
      fmt_s = imformats(extension);
      tf = feval(fmt_s.isa, filename);

      % If we retrieve something, extract the actual full extension
      if (tf)
        format = fmt_s.ext{1};
      end
    end

    % In case we do not have anything, go for the slow version
    if (isempty(format))
      infos = imfinfo(fname);

    % Otherwise, we can directly call the adequat parsing function
    else
      infos = feval(fmt_s.info, filename);
    end

    % Finally we can get the number of frames and the image size
    nframes = length(infos);
    ssize = [infos(1).Height infos(1).Width];
  end

  return;
end
