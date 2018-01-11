function [p, fval] = myfit(varargin)

  fitting = parse_inputs(varargin{:});

  ndecimals = -min(floor(real(log10(fitting.tolerance/10))), 0);
  if (~isfinite(ndecimals))
    ndecimals = 0;
  end

  nparams = length(fitting.init_pos);
  p0 = fitting.init_pos .* (1 + randn([nparams,1])*fitting.init_noise);

  scaling_functions = cell(nparams, 1);
  bounded = (isfinite(fitting.bounds));

  x0 = p0 ./ fitting.init_pos;
  for i=1:nparams
    if (all(bounded(i,:)))
      scaling_functions{i} = @ul_finite;
      fitting.bounds(i,2) = diff(fitting.bounds(i,:));

      x0(i) = (x0(i) - fitting.bounds(i,1)) / fitting.bounds(i,2);
    elseif (bounded(i,1))
      scaling_functions{i} = @l_finite;

      x0(i) = (x0(i) - fitting.bounds(i,1)) / (x0(i) - fitting.bounds(i,1) + 1);
    elseif (bounded(i,2))
      scaling_functions{i} = @u_finite;

      x0(i) = 1 / (fitting.bounds(i,2) - x0(i) - 1);
    else
      scaling_functions{i} = @none_finite;

      x0(i) = (x0(i) - 2 + sqrt(x0(i)^2 + 4) ) / (2* x0(i));
    end
  end

  opt = cmaes('defaults');
  opt.MaxFunEvals = fitting.max_iter;
  opt.TolFun = fitting.tolerance;
  opt.TolX = fitting.tolerance/10;
  opt.SaveFilename = '';
  opt.SaveVariables = 'off';
  opt.EvalParallel = 'yes';
  opt.LogPlot = 0;
  opt.LogFilenamePrefix = fitting.logging_name;
  opt.StopOnWarnings = false;
  opt.WarnOnEqualFunctionValues = false;
  opt.PopSize = min(10*numel(p0), fitting.max_population);
  opt.DiagonalOnly = (numel(p0) > 100);

  [x, fval, cmaes_count, stopflag, out] = cmaes(@error_function, x0(:), fitting.step_size, opt);
  x = out.solutions.bestever.x;
  fval = out.solutions.bestever.f;

  options = optimset('MaxFunEvals', max(fitting.max_iter - cmaes_count,0), ...
                     'MaxIter', max(fitting.max_iter - cmaes_count,0), ...
                     'Display', 'off', ...
                     'TolFun', fitting.tolerance, ...
                     'TolX', fitting.tolerance/10);

  if (~isempty(fitting.logging_name))
    options.OutputFcn = @log_simplex;
    fid = fopen([fitting.logging_name 'evol.dat'], 'a');
  end

  display('Local optimization using the Simplex method');
  [x, fval_opt, stopflag_opt, out_opt] = fminsearch(@error_function, x(:), options);

  if (~isempty(fitting.logging_name))
    fclose(fid);
  end

  disp(['Simplex: ' num2str(fval) ' -> ' num2str(fval_opt)]);
  fval = fval_opt;

  x = roundn(x, -ndecimals);
  p = convert_parameters(x, fitting.bounds, fitting.init_pos, scaling_functions);

  return;

  function stop = log_simplex(x, optimValues, state)

    stop = false;

    if ((mod(optimValues.iteration, 10) == 1 && strncmp(state, 'iter', 4)) || strncmp(state, 'done', 4))
      print_str = [' %.' num2str(ndecimals) 'f'];

      if (uuid(1) == 'C')
        fprintf(fid, [uuid '1 %ld  %e 0 0'], optimValues.iteration+cmaes_count, optimValues.fval);
        fprintf(fid, print_str, x);
        fprintf(fid, '\n');
      else
        fprintf(fid, [uuid '%ld : %ld |'], optimValues.iteration+cmaes_count, optimValues.fval);
        fprintf(fid, print_str, x);
        fprintf(fid, '\n');
      end

      disp([num2str(optimValues.iteration+cmaes_count) ' : ' num2str(optimValues.fval), ' | ' num2str(x(:).')]);
    end

    return;
  end

  function err_all = error_function(varargin)

    p_all = roundn(varargin{1}, -ndecimals);

    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams;
    end
    err_all = NaN(1, nevals);

    outs = any(p_all<0 | p_all>1, 1);

    p_scaled = convert_parameters(p_all, fitting.bounds, fitting.init_pos, scaling_functions);

    for i = 1:nevals
      if (~outs(i))
        err_all(i) = fitting.error_function(p_scaled(:,i), fitting.error_parameters);
      end
    end

    err_all(isnan(err_all)) = Inf;

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = fitting.max_error*sign(err_all(1));
    end

    return;
  end
end

function x = convert_parameters(x, bounds, scaling, funcs)

  for i=1:size(x,1)
    x(i,:) = scaling(i)*funcs{i}(x(i,:), bounds(i,:));
  end

  return;
end

function p = ul_finite(x, bounds)

  p = x*bounds(2) + bounds(1);

  return;
end

function p = u_finite(x, bounds)

  p = bounds(2) - (1 - x) ./ x;

  return;
end

function p = l_finite(x, bounds)

  p = x ./ (1 - x) + bounds(1);

  return;
end

function p = none_finite(x, bounds)

  p = x ./ (1 - x) - (1 - x) ./ x;

  return;
end

function fitting = parse_inputs(varargin)

  fitting = get_struct('fitting');

  skip = false;
  for i=1:length(varargin)
    if (skip)
      skip = false;
    else
      if (isstruct(varargin{i}))
        fitting = varargin{i};
      elseif (isnumeric(varargin{i}))
        if (isempty(fitting.init_pos))
          fitting.init_pos = varargin{i};
        else
          fitting.bounds = varargin{i};
        end
      elseif (ischar(varargin{i}))
        if (i == length(varargin))
          error('Properties pairs must come in PAIRS.');
        else
          % If the parameter exists in opts we simply assign it the
          % provided value
          if (isfield(fitting, varargin{i}))
            fitting.(varargin{i}) = varargin{i+1};

          % Or ignore it !
          else
            warning(['Property ''' varargin{i} ''' does not exist. Ignoring']);
          end

          skip = true;
        end
      elseif (isa(varargin{i}, 'function_handle'))
        fitting.error_function = varargin{i};
      else
        warning(['Unknown type of input ''' class(varargin{i}) '''. Ignoring']);
      end
    end
  end

  fitting.init_pos = fitting.init_pos(:);

  if (isempty(fitting.error_function))
    error('No error function specified.');
  elseif (isempty(fitting.init_pos))
    error('No initial guess for the parameters specified.');
  elseif (isempty(fitting.bounds))
    fitting.bounds = Inf(length(fitting.init_pos), 2);
    fitting.bounds(:,1) = -fitting.bounds(:,1);
  end

  if (any(fitting.bounds(:,1) >= fitting.bounds(:,2)))
    error('All defined intervals for the values must be valid.');
  end

  return;
end
