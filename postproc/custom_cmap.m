


%% Define custom colormap

D=[1 1 1;
   0 0 1;
   0 1 0;
   1 1 0;
   1 0 0;];
F=[0 0.25 0.5 0.75 1];
G=linspace(0,1,256);
cmap=interp1(F,D,G);


wbgyr = cmap;


D=[0 0 1;
   1 1 1;
   1 0 0;];
F=[0 0.5 1];
G=linspace(0,1,256);
bwr=interp1(F,D,G);


