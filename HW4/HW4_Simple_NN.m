vars = load('simple_nn_vars.mat');
n_epoch = 2000;
rate = [0.04, 0.2];
test_intv = 100;

output_nn = simple_nn(vars, n_epoch, rate, test_intv);
%

function output = simple_nn(vars, n_epoch, rate, test_intv)
%SIMPLE_NN Performs simple neural network estimation of linear model
%    INPUTS
%     vars: Structure containing relevant variables
%  n_epoch: number of epochs to compute
%     rate: learning rate
%test_intv: number of training epochs to pass between testing validation
%   OUTPUTS
%   output: structure containing error and classification results

% Extract variables
system_len  = vars.system_len;
test_in     = vars.test_in;
test_out    = vars.test_out;
train_in    = vars.train_in;
train_out   = vars.train_out;

% Initialize random seed
rng(0);

% Initialize synapses
syn0 = 2.*rand(system_len, 4) - 1;
syn1 = 2.*rand(4, 1) - 1;

% Initialize error arrays and test iterator
l2_err_train = zeros(n_epoch, 1);
l2_err_test = zeros(floor(n_epoch/test_intv) + 1, 2);
ts_i = 1;

% Run neural network training/testing for given number of epochs
for i = 1:n_epoch
    % Run training forward propagation
    [l0, l1, l2, l2_err] = simple_nn_fwd(train_in, train_out, syn0, syn1);
    
    % Compute mean squared testing error for final layer
    l2_err_train(i) = mean(l2_err.^2);
    
    % Regularly get nnet accuracy for test data
    if (mod(i, test_intv) == 0) || (i == 1)
        % Run testing forward propagation
        [~, ~, ~, l2_ts_err] = simple_nn_fwd(test_in, test_out, syn0, syn1);
        
        % Compute mean squared error for final layer and iterate index
        l2_err_test(ts_i,:) = [i, mean(l2_ts_err.^2)];
        ts_i = ts_i + 1;
    end
    
    % Perform backpropagation
    [syn0, syn1] = backprop_simple(l0, l1, l2, l2_err, syn0, syn1, rate);
end

% Build output structure
output.syn0 = syn0;
output.syn1 = syn1;
output.l2_err_train = l2_err_train;
output.l2_err_test = l2_err_test;

end

%------------------------------------------------------------------------------%
%---------------------------Add Helper Functions Here--------------------------%
%------------------------------------------------------------------------------%
function [l0, l1, l2, l2_err] = simple_nn_fwd(train_in, train_out, syn0, syn1)
% Initial layer
l0 = train_in;
% Middle layer
suml0 = l0*syn0;
l1 = sigmoid(suml0,0);
% Output layer
suml1 = l1*syn1;
l2 = sigmoid(suml1,0);
l2_err = train_out - l2;
end

function [syn0, syn1] = backprop_simple(l0, l1, l2, l2_err, syn0, syn1, rate)
% Output layer
deltal2 = l2_err.*sigmoid(l2,1);
l1_err = deltal2*syn1';
% Middle layer
deltal1 = l1_err.*sigmoid(l1,1);
syn0 = syn0 + rate(1)*l0'*deltal1;
syn1 = syn1 + rate(2)*l1'*deltal2;
end

function output = sigmoid(x, dodv)
sig = 1./(1+exp(-x));
if dodv == 0
    output = sig;
elseif dodv == 1
    output = (1-x).*x;
end
end