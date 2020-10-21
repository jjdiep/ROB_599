%% Script to run separate HW4 parts

%% Run P2
vars = load('simple_nn_vars.mat');
n_epoch = 2000;
rate = [0.04, 0.2];
test_intv = 100;

output_nn = simple_nn(vars, n_epoch, rate, test_intv);

%% Run P3
vars = load('simple_cnn_vars.mat');
n_epoch = 700;
rate = 0.01;
p = 2;
test_intv = 50;

output_cnn = simple_cnn(vars, n_epoch, rate, p, test_intv);