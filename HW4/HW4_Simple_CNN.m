vars = load('simple_cnn_vars.mat');
n_epoch = 700;
rate = 0.01;
p = 2;
test_intv = 50;

output_cnn = simple_cnn(vars, n_epoch, rate, p, test_intv);
%
function output = simple_cnn(vars, n_epoch, rate, p, test_intv)
%SIMPLE_CNN Performs simple CNN classification of two types of symbols
%    INPUTS
%     vars: Structure containing relevant variables
%  n_epoch: number of epochs to compute
%     rate: learning rate
%        p: pooling size/step
%test_intv: number of training epochs to pass between testing validation
%   OUTPUTS
%   output: structure containing error and classification results

% Extract inputs from struct
train_cell      = vars.train_cell;
train_labels    = vars.train_labels;
test_cell       = vars.test_cell;
test_labels     = vars.test_labels;
k_cell          = vars.k_cell;
weights         = vars.weights;

% Build arrays to store errors and number of correct classifications
train_epoch_err = zeros(n_epoch, 1);
tr_num_corr_arr = zeros(n_epoch, 1);
test_epoch_err  = zeros(floor(n_epoch/test_intv) + 1, 2);
ts_num_corr_arr = zeros(floor(n_epoch/test_intv) + 1, 2);
ts_i = 1;   % Test iterator

% Because we are only training final weights, pre-compute FCNs
train_fconv = cell(size(train_cell));
test_fconv = cell(size(test_cell));
for i = 1:numel(train_cell)
    train_fconv{i} = cnn_fwd(train_cell{i}, k_cell, p);
end
for i = 1:numel(test_cell)
    test_fconv{i} = cnn_fwd(test_cell{i}, k_cell, p);
end


% Run cnn training/testing
for i = 1:n_epoch
    % Initialize training error and correct classification vars for iteration
    train_err = zeros(numel(train_cell), 1);
    tr_num_corr = 0;
    
    % Train cnn for each training image
    for j = 1:numel(train_cell)
        % Pull training fully connected layer
        fconv = train_fconv{j};
        
        % Apply backpropagation, compute mean-squared error
        [weights, tn_err] = backprop_cnn(fconv,weights,train_labels(j,:),rate);
        train_err(j) = sum(tn_err.^2);
        
        % Check training classifications
        tr_guess = fconv*weights;% complete ...
        tr_num_corr = tr_num_corr + eval_class(tr_guess, train_labels(j,:));
    end
    
    % Regularly compute testing performance
    if (mod(i, test_intv) == 0) || (i == 1)
        % Initialize testing error and correct classification vars
        test_err = zeros(numel(test_cell), 1);
        ts_num_corr = 0;
        
        % Evaluate accuracy with testing set (no backprop)
        for j = 1:numel(test_cell)
            % Pull testing fully connected layer
            fconv = test_fconv{j};
            
            % Compute estimate for testing image and compute mean-squared error
            ts_guess = fconv*weights; % complete ...
            ts_err =  test_labels(j,:)-ts_guess; % complete ...
            test_err(j) = sum(ts_err.^2);
            
            % Check testing classifications
            ts_num_corr = ts_num_corr + eval_class(ts_guess, test_labels(j,:));
        end
        
        % Store testing errors/classification accuracy and epochs
        test_epoch_err(ts_i,:) = [i, mean(test_err)];
        ts_num_corr_arr(ts_i,:) = [i, ts_num_corr];
        ts_i = ts_i + 1;    % Iterate testing evaluation index
    end
    
    % Store training errors/classification accuracy
    train_epoch_err(i) = mean(train_err);
    tr_num_corr_arr(i) = tr_num_corr;
end

% Build output structure
output.train_epoch_err = train_epoch_err;
output.test_epoch_err  = test_epoch_err;
output.tr_num_corr_arr = tr_num_corr_arr;
output.ts_num_corr_arr = ts_num_corr_arr;
output.weights = weights;

% Plots
figure()
hold on
plot(1:n_epoch, train_epoch_err, 'k.-')
plot(test_epoch_err(:,1), test_epoch_err(:,2), 'bo-')
hold off
xlabel('Epoch')
ylabel('Mean Sq. Error')
legend('Training Error', 'Testing Error', 'Location', 'northeast')

tr_n = numel(train_cell);
ts_n = numel(test_cell);
figure()
hold on
plot(1:n_epoch, 100*(tr_n - tr_num_corr_arr)./tr_n, 'k.-')
plot(ts_num_corr_arr(:,1), 100*(ts_n - ts_num_corr_arr(:,2))./ts_n, 'bo-')
hold off
ylim_cur = ylim;
ylim([0, ylim_cur(2)])
xlabel('Epoch')
ylabel('Classification Loss (%)')
legend('Training Loss', 'Testing Loss', 'Location', 'northeast')

end

% Change this for CNN
function fconv = cnn_fwd(mod_cell, k_cell, p)
for i = 1:numel(k_cell)
    conv_layer = conv2(mod_cell,k_cell{i},'same');
    relu_layer = relu(conv_layer);
    pool_layer = pool(relu_layer,p);
    
    conv_layer2 = conv2(pool_layer,k_cell{i},'same');
    relu_layer2 = relu(conv_layer2);
    pool_layer2 = pool(relu_layer2,p);
    
    conv_layer3 = conv2(pool_layer2,k_cell{i},'same');
    relu_layer3 = relu(conv_layer3);
    pool_layer3 = pool(relu_layer3,p);
    
    norm_factor = sum(sum(pool_layer3));
    if norm_factor == 0
        norm_layer = pool_layer3;
    else
        norm_layer = pool_layer3/norm_factor;
    end
    foldconv(4*(i-1)+1:4*(i-1)+4,1) = [norm_layer(:,1); norm_layer(:,2)];
end
fconv = foldconv';
% for j = 1:1:numel(foldconv)
%     fconv(1,j) = foldconv(j);
% end
% fconv = ([foldconv(:,1); foldconv(:,2)])';
end

function [weights, tn_err] = backprop_cnn(fconv,weights,train_labels,rate)
result = fconv*weights;
err = train_labels - result;
rel = rate*err'*fconv;
weights = weights + rel';
tn_err = train_labels - fconv*weights;
end

function accuracy = eval_class(guess,label)
if guess(1) > guess(2)
    myguess = [1 0];
else
    myguess = [0 1];
end
if myguess == label
    accuracy = 1;
else
    accuracy = 0;
end
end

function output = relu(x)
output = max(0,x);
end

function p_layer = pool(arr,p)
if mod(length(arr),p) ~= 0
    modcalc = mod(length(arr),p);
    addzeros = p - modcalc;
    arr = [arr zeros(length(arr),addzeros); zeros(addzeros,length(arr)+addzeros)];
end
for i = 1:1:length(arr)/p
    for j = 1:1:length(arr)/p
        newarr = arr((p*(i-1)+1):p*i,(p*(j-1)+1):p*j);
        out = max(newarr);
        p_layer(i,j) = max(out);
    end
end
end