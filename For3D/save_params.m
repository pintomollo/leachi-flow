function save_params(name_parameters, params, filename, write_mode)

%% function save_params(name_parameters, params, filename, write_mode)
%
% Save parameters name and value in a file named params.txt by default
% at the begining of each For3D run.
%
% -- inputs --
% * name_parameters: description of each parameter
% * params: value of each parameter
% * filename: name of the saved file (default, params.txt)
% * write_mode: write or append (default)
%
% -- output --
% The list of parameters description and value is automatically saved
% in the current folder.
%
% (c) Arnauld SERGE, 2011
%
% see also rendering_3D, FOR3D_dialog_box


if nargin<3, filename = 'params.txt'; end
if nargin<4, write_mode = 'at'; end % wt??

%   while fid==-1,
fid = fopen(filename, write_mode); %, 'wt','native');
%       pause (.01)
%   end

if any(write_mode ~= 'a')
    fprintf(fid, '%s @ %s ', filename, date) ;
    clk = clock() ;
    fprintf(fid, '%.2dh%.2dm%.2ds\n\n', clk(4), clk(5), round(clk(6)));
end

for i=1:length(params)
    if iscell(params{i})
        params{i} = cell2mat(params{i});
    end
    if ischar(params{i})
        template = '%s : %s\n';
    else
        template = '%s : %g\n';
    end
    fprintf(fid, template, name_parameters{i}, params{i});
end

fclose(fid);

%%%