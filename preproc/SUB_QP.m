clear all;

d = pwd;

MAKE_RP;

cd('/u/scratch/s/sgess');
stat = mkdir(rpinput_output_name);
cd(rpinput_output_name);
copyfile('/u/home/mori/sgess/executables/QuickPIC/qpic.e.twiss.0907','qpic.e');
copyfile([d '/rpinputs/' date_dir 'rpinput_' rpinput_output_name],'rpinput');
copyfile([d '/params/' date_dir 'param_' rpinput_output_name '.mat'],'param.mat');
copyfile([d '/commands/' date_dir 'qpic.e.cmd_' rpinput_output_name],'qpic.e.cmd');

%unix('qsub qpic.e.cmd','-echo');

cd(d);
exit;
