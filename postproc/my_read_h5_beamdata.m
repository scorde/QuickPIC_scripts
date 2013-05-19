function beam = my_read_h5_beamdata(filename)

char_slash = findstr(filename, '/');
char_dot = findstr(filename, '.');

dataname = filename(char_slash(end):char_dot(end)-1);
dataname( dataname == '_' ) = '-'; % for some reason QuickPIC has slightly different filename than dataname

%fileinfo = h5info(filename);
% Read a subset of the data using info structure

%p1 = h5read(filename, '/p1');
%p2 = h5read(filename, '/p2');
%p3 = h5read(filename, '/p3');
%x1 = h5read(filename, '/x1');
%x2 = h5read(filename, '/x2');
%x3 = hf5read(filename, '/x3');

% OLD LAPTOP MATLAB VERSION
p1 = hdf5read(filename, '/p1');
p2 = hdf5read(filename, '/p2');
p3 = hdf5read(filename, '/p3');
x1 = hdf5read(filename, '/x1');
x2 = hdf5read(filename, '/x2');
x3 = hdf5read(filename, '/x3');

beam = [x1 x2 x3 p1 p2 p3];

%hist(x3, 51)
%stop

return;

