function [choc,channel,plateau] = ind2tab(index)

% Get choc, channel and plateau number from the index of spectrum

load('database_index.mat');

choc = database_index.index_choc(index);
channel = database_index.index_channel(index);
plateau = database_index.index_plateau(index);

end
