function remove_3D_corner(filename)

%% function remove_3D_corner(filename)
%
% Remove top/front corner from a 3D stack (according to the default 3D orientation):
% data above Nx/2, Ny/2 & Nz/2 are set to 0, allowing to view inside the volume.
%
% -- input --
% filename: image stack name, expected to be already filtered by a
% previous run of For3D, hence named with the extention '_filt.tif'
%
% -- output --
% The resulting stack is automatically saved, with the extension '_no_corner.tif'
% (or '_no_corner_filt.tif' if the input is already filtered).
%
% (c) AS 2014
%
% see also plot_3D_corner, rendering_3D


if nargin < 1, filename = dir('*_filt.tif'); end

files = dir(filename);
Nf = length(files);

for nf = 1:Nf
    file = files(nf).name;
    
    if strfind(file, '_filt.tif'), no_corner_file = [file(1:end-9) '_no_corner_filt.tif'];
    else no_corner_file = [file(1:end-4) '_no_corner.tif'];
    end
    
    im = imread(file, 'Index', 1);
    Nx = size(im, 1);
    Ny = size(im, 2);
    Nz = length(imfinfo(file));
    
    fprintf('loading and removing corner from %d images in %s    ', Nz, file)
    for nz = 1:Nz
        im = imread(file, 'Index', nz);
        
        if nz >= ceil(Nz/2)
            im(1:floor(Nx/2), 1:floor(Ny/2), :, :) = 0;
        end
        
        if (nz==1), imwrite(uint8(im), no_corner_file, 'tiff', 'Compression', 'none')
        else imwrite(uint8(im), no_corner_file, 'tiff', 'Compression', 'none', 'WriteMode', 'append')
        end
        
        fprintf('\b\b\b%3d', nz)
    end
    fprintf('\n')
end

%%%