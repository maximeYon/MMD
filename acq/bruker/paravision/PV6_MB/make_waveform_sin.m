% Sin gradient waveform
%
% Output in the format for Matthew Budde's code for Bruker pv6.0.1 

% Define path for output folders
%[run_path,run_name,run_ext] = fileparts(run_fn);
run_path = cd;
out_path = fullfile(run_path,'waveforms');


% Parameters for calculating b-values
gmax = 100/100*0.6; % Max gradient [3 T/m]
tau = 10e-3; % Waveform duration [10 ms]

% Define timing parameters relative to a total echo time = 1
np = 200; % Number of time intervals in waveform [1000]
epsilon_up = .5; % Quarter-sine ramp up [0.5]
plateau = 0; % Plateau [0]
epsilon_down = .5; % Quarter-sine ramp down [0.5]
out_fn = fullfile(out_path,'sin_200');

%------------------------

taured = 1;
dt = taured/np;
t = taured*linspace(0,1,np);

% Quarter-sine ramp up
np_epsilon_up = round(epsilon_up/dt);
t_epsilon_up = pi/2*linspace(0,1,np_epsilon_up)';
g_up = sin(t_epsilon_up);
%figure(1), clf, plot(t_epsilon_up,g_up,'-'), return

% Quarter-sine ramp down
np_epsilon_down = round(epsilon_down/dt);
t_epsilon_down = pi/2*linspace(1,2,np_epsilon_down)';
g_down = sin(t_epsilon_down);
%figure(1), clf, plot(t_epsilon_down,g_down,'-'), return

% Plateau
np_plateau = round(plateau/dt);

ga = [g_up; ones(np_plateau,1); g_down];
ga = mean([ga flipud(ga)],2); % Removes asymmetry
ga([1 end]) = 0;
%figure(1), clf, plot(t,ga,'-'), return

gx = 0*ga;
gy = 0*ga;
gz = ga;

figure(1), clf
plot(t,gx,'r-',t,gy,'g-',t,gz,'b-')
ylabel('g / gmax')
xlabel('t / \tau')

% Save waveforms in Bruker format
[out_path,out_name,out_ext] = fileparts(out_fn);
if ~isdir(out_path)
    mkdir(out_path)
end

title_str = ['Quarter-sine ramp up and down'];

fid = fopen(out_fn,'w');

text.header = {
{['##TITLE= ' title_str]};
{'##JCAMP-DX= 5.00 Bruker JCAMP library'}
{'##DATA TYPE= Shape Data'}
{'##ORIGIN= Bruker Analytik GmbH'}
{'##OWNER= <nmrsu>'}
{['##DATE= ' datestr(now,'yyyy-mm-dd')]}
{['##TIME= ' datestr(now,'HH:MM:SS')]}
{'##MINX= 0'}
{'##MAXX= 1'}
{'##MINY= 0'}
{'##MAXY= 1'}
{'##$SHAPE_EXMODE= Gradient'}
{'##$SHAPE_TOTROT= 0'}
{'##$SHAPE_BWFAC= 0'}
{'##$SHAPE_INTEGFAC= 0'}
{'##$SHAPE_MODE= 0'}};

[nlines, ~] = size(text.header);

for nline = 1:nlines
    fprintf(fid,'%s\n',text.header{nline}{1});
end

fprintf(fid,'%s\n',['##BFACTOR=']);
fprintf(fid,'%s\n',['##DURATION=' num2str(tau)]);
fprintf(fid,'%s\n',['##DIRECTIONVEC= 0 0 1']);
fprintf(fid,'%s\n',['##NPOINTS=' num2str(numel(gx))]);
fprintf(fid,'%s\n',['##XYDATA= (T X Y Z)']);

formatspec = '%8.6f %8.6f %8.6f\r\n';

out_mat = [gx gy gz];
fprintf(fid, formatspec, out_mat');

fprintf(fid,'%s\n',['##END']);
fclose(fid);
