function [axis_id] = plot_time_mf (nuc_name, xmf, time, varargin)
%--------------------------------------------------------------------------
%[] = plot_time_mf (nuc_name, xmf, time, ...)
% Plots mass fractions vs. time for provided data. 
% Inputs>  nuc_name: array of species names 
%          xmf: time dependent array of mass fractions 
%          time: array of times
% Options: SeparateLegend: draw separate legends for each lineset
%          LineStyle: line style of mass fraction lines
%          LineWidth: width of mass fraction lines
%          TimeFromEnd: plot time relative to end of calculation
% Outputs< axis_id: handle of current axis
%--------------------------------------------------------------------------

% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('nuc_name',@iscellstr);
  p.addRequired('xmf',@(x)validateattributes(x, {'numeric'}, {'2d', 'real', ...
    '>=', 0, '<=', 1.1}));
  p.addRequired('time',@(x)validateattributes(x, {'numeric'}, {'vector', ...
    'real', '>=', 0}));

% Define optional inputs
  p.addOptional('SeparateLegend', false, @(x)validateattributes(x, {'logical'},{'scalar'}));
  p.addOptional(     'LineStyle',    '', @(x)validateattributes(x, {'char'},{'row'}));
  p.addOptional(     'LineWidth',     2, @(x)validateattributes(x, {'numeric'}, ...
    {'scalar', 'real', 'positive', '>=', 0.1, '<=', 10}))
  p.addOptional(   'TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{'scalar'}));

% Parse and validate all input arguments.
  p.parse(nuc_name, xmf, time, varargin{:});

% Assign input arguments
  nuc_name        = p.Results.nuc_name;
  xmf             = p.Results.xmf;
  time            = p.Results.time;
  separate_legend = p.Results.SeparateLegend;
  linestyle       = p.Results.LineStyle;
  linewidth       = p.Results.LineWidth;
  time_from_end   = p.Results.TimeFromEnd;

%Set default line colors and types, if not over
  set(0,'DefaultAxesColorOrder',[[0 0 0];[1 0 0];[1 .5 0];[1 1 0];[.5 1 0];[0 1 0];[0 1 .5];[0 1 1];[0 .5 1];[0 0 1];[.5 0 1];[1 0 1];[1 0 .5];[.5 .5 .5]]);
  set(0,'DefaultAxesLineStyleOrder','-|-.|--|:')

% Test dynamic range in time
  ntime=size(time,2);
  if( time_from_end == true);
    time_stop=time(1);
    time_start=time(ntime-1);
    time_range=time_start/time_stop;
     time_label='Time from Completion (s)';
    time_direction='reverse';
 else
    time_start=time(2);
    time_stop=time(ntime);
    time_range=time_stop/time_start;
    time_label='Time(s)';
    time_direction='normal';
  end
 
% For small temporal range, use linear time axis 
  if(time_range > 20) 
    plot_id=loglog(time,xmf,linestyle,'LineWidth',linewidth);

% For large temporal range, use log time axis 
  else      
    plot_id=semilogy(time,xmf,linestyle,'LineWidth',linewidth);
    
  end

% Label and limit Plot  
  axis_id = gca; 
  ylabel('Mass fraction');
  xlabel(time_label);
  set(axis_id,'XDir',time_direction);
  set(axis_id,'XLim',[time_start time_stop]);
  
% Build Legend
  legend_array=cell(nuc_name);
  [legend_id,hobjs,hlinpatches,legtxt] = legend(plot_id,char(legend_array),'Location','SouthEast');
  if(separate_legend)
    hleg_copy = copyobj(legend_id,gcf);
    delete(legend_id);
  end
  
end

