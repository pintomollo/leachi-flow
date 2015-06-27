function pipeline(varargin)

  disp('----B. leachi pipeline----')
  disp('Parsing inputs...')
  [myrecording, mysimulation, opts] = parse_inputs(varargin{:});

  do_simulate = false;
  if (isempty(myrecording.channels))
    disp('Simulating flow...')
    [fname, mysimulation] = simulate_flow(mysimulation);
    opts.simulation = mysimulation;

    myrecording.channels(1).fname = fname;
    [junk, exp_name, junk] = fileparts(fname);
    [junk, exp_name, junk] = fileparts(exp_name);
    myrecording.experiment = exp_name;

    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
    do_simulate = true;
  end

  exp_name = myrecording.experiment;

  do_process = false;
  for i=1:length(myrecording.channels)
    if (isempty(myrecording.channels(i).file))
      disp('Converting recording...')

      if (do_simulate)
        myrecording.channels(i).fname = convert_movie(myrecording.channels(i).fname, false);
      else
        myrecording.channels(i).fname = convert_movie(myrecording.channels(i).fname);
      end
      do_process = true;

      if (isempty(exp_name))
        exp_name = myrecording.channels(i).fname;
        [junk, exp_name, junk] = fileparts(exp_name);
        [junk, exp_name, junk] = fileparts(exp_name);
        exp_name = regexprep(exp_name, ' ', '');
        myrecording.experiment = exp_name;
      end

      save([myrecording.experiment '.mat'], 'myrecording', 'opts');
    end
  end

  if (do_process)
    disp('Processing recording...')
    [myrecording, opts] = preprocess_movie(myrecording, opts);
    save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  end

  disp('Analyzing flow...')
  myrecording = leachi_flow(myrecording, opts);
  save([myrecording.experiment '.mat'], 'myrecording', 'opts');

  disp('DONE !')

  return;
end

function [myrecording, mysimulation, opts] = parse_inputs(varargin)

  myrecording = get_struct('myrecording');
  mysimulation = get_struct('simulation');
  opts = get_struct('options');

  conf_rec = '';
  conf_sim = '';
  conf_opt = '';

  % First run: get structures, load config files & remove them
  for i=length(varargin):-1:1
    if (isstruct(varargin{i}))
      if (isfield(varargin{i}, 'experiment'))
        myrecording = varargin{i};
        varargin(i) = [];
      elseif (isfield(varargin{i}, 'init_simulation'))
        mysimulation = varargin{i};
        varargin(i) = [];
      elseif (isfield(varargin{i}, 'config_files'))
        opts = varargin{i};
        varargin(i) = [];
      end
    elseif (ischar(varargin{i}))
      switch varargin{i}
        case 'config_recording'
          conf_rec = varargin{i+1};
          varargin(i:i+1) = [];
        case 'config_simulation'
          conf_sim = varargin{i+1};
          varargin(i:i+1) = [];
        case 'config_options'
          conf_opt = varargin{i+1};
          varargin(i:i+1) = [];
      end
    end
  end

  % Retrieve the latest version of the structures
  myrecording = update_structure(myrecording, 'myrecording');
  mysimulation = update_structure(mysimulation, 'simulation');
  opts = update_structure(opts, 'options');

  % Load the configuration files
  if (~isempty(conf_rec))
    myrecording = load_parameters(myrecording, conf_rec);
  end
  if (~isempty(conf_sim))
    mysimulation = load_parameters(mysimulation, conf_sim);
  end
  if (~isempty(conf_opt))
    opts = load_parameters(opts, conf_opt);
  end

  % Now we check that the parameters were provided in pairs
  npairs = length(varargin) / 2;
  if (npairs ~= floor(npairs))
    error('Properties pairs must come in PAIRS.');
  end

  % Second run: loop over the pairs of parameters and assign them
  for i = 1:npairs

    % Get the pair
    field = varargin{2*i - 1};
    value = varargin{2*i};

    % Try to assign them to the various structures
    [opts, is_valid] = assign_field(opts, field, value);
    if (~is_valid)

      [myrecording, is_valid] = assign_field(myrecording, field, value);

      if (~is_valid)
        [mysimulation, is_valid] = assign_field(mysimulation, field, value);
      end
    end

    % Otherwise we notify its inexistence
    if (~is_valid)
        warning(['Property ''' field ''' does not exist. Ignoring']);
    end
  end

  return;
end
