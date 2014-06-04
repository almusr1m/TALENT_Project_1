function [] = compare_time_mf( filenames, varargin )
%--------------------------------------------------------------------------
%[] = compare_time_mf ( filenames, mass_limit, max_sort, ...)
% Plots mass fractions vs. time for each filename. 
% Inputs>  filenames: a cell array of file names {'file1'; 'file2'; ...}
%          mass_limit: a limiting mass fraction to include in plot
%          max_sort: when true, species are ordered by maximium mass 
%            fraction; else, lines are ordered as in data file. 
% Outputs< None
%--------------------------------------------------------------------------

% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('filenames',@iscellstr);

% Define optional inputs
  p.addOptional(  'MassLimit', 1e-25, @(x)validateattributes(x, {'numeric'}, ...
                {'scalar', 'real', 'positive', '>=', 1e-30, '<=', 1e-1}));
  p.addOptional(    'MaxSort',  true, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional('TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional(     'ZLimit', [0 0], @(x)validateattributes(x, {'numeric'}, ...
    {'vector', 'integer'}));

% Parse and validate all input arguments.
  p.parse(filenames, varargin{:});

% Assign input arguments
  filenames     = p.Results.filenames;
  mass_limit    = p.Results.MassLimit;
  max_sort      = p.Results.MaxSort;
  time_from_end = p.Results.TimeFromEnd;
  z_limit       = p.Results.ZLimit;
  
% Define manual linestyles  
  linestyles = {'-';'--';':';'-.'}
  linestyles(3)
  figure;

% Loop over files  
  num_files = size(filenames,1);
  for ifile= num_files:-1:1;
    filename = char(filenames(ifile,:));
  
% Choose file type
    is_ev_file=strfind(filename,'ev');
    if(isempty(is_ev_file));
        
% Read TS file
      [zz, aa, xmf, time, ~, ~, ~, ~, ~ , ~] = read_ts_file(filename);

% Build Isotope symbols
      [ nuc_name ] = build_isotope_symbol ( zz,aa );

% Limit abundances plotted by element (zmin < Z <zmax)
    if (z_limit(1)~=0)
      zmin = z_limit(1);
    else
      zmin = min(zz);
    end
    if (z_limit(2)~=0)
      zmax = z_limit(2);
    else
      zmax = max(zz);
    end
    limit_z = find(zz >= zmin & zz<=zmax)
    xmf=xmf(limit_z,:);
    nuc_name=nuc_name(limit_z);

    else  
% Alternately, the ev_file may be used
      [ nuc_name, xmf, time, ~, ~, ~, ~] = read_ev_file (filename );

    end

% Limit abundances plotted by maximum mass
    xmf_max= max(xmf,[],2);
    limit_xmf = find(xmf_max > mass_limit);
    xmf=xmf(limit_xmf,:);
    nuc_name=nuc_name(limit_xmf);
  
% Sort isotopes in descending order of maximum mass fraction
    if (max_sort == true);
      xmf_max= max(xmf,[],2);
      [xmf_sort,sort_order]  = sort(xmf_max,'descend');
      nuc_name=nuc_name(sort_order);
      xmf=xmf(sort_order,:);
    end
    
% Plot Mass Fraction
    if (ifile > 1) 
      separate_legend = true
    else
      separate_legend = false
    end
    linestyle=char(linestyles(ifile))
    size(linestyle)
    mf_axis_id=plot_time_mf (nuc_name, xmf, time, 'LineWidth', 2, ...
      'LineStyle', linestyle, 'SeparateLegend', separate_legend,...
      'TimeFromEnd',time_from_end); 
    set(mf_axis_id,'Ylim',[mass_limit 2])
    hold on
  end
  hold off
end

