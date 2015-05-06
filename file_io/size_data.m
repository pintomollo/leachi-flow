function [nframes, ssize] = size_data(fname)
% SIZE_DATA extracts the number of frames as well as the size of a frame from a TIFF file.
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
      fname = fname.fname;
      fid = fopen(fname, 'r');
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

  % Now we can work !
  else

    % Get the full filename just in case
    filename = fopen(fid);
    fclose(fid);

    % Use the Tiff library libtiff using the gateaway provided by Matlab
    try
      tif = Tiff(filename, 'r');

    % Just in case we try to open something else
    catch
      warning('Tracking:SizeData', '%s is not a compatible TIFF file.', filename);
      return;
    end

    % Count the number of frames by going through them
    nframes = 1;
    while (~tif.lastDirectory())
      nframes = nframes + 1;
      tif.nextDirectory();
    end

    % Extract the size of the image
    ssize = [tif.getTag('ImageLength') tif.getTag('ImageWidth')];

    % And close it
    tif.close();
  end

  return;
end
