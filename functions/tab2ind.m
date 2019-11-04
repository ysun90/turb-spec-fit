function index = tab2ind(choc,channel,plateau)

% Get index of spcetrum from choc, channel and plateau

load('database_index.mat');

index = find(database_index.index_choc==choc & ...
             database_index.index_channel==channel & ...
             database_index.index_plateau==plateau);
    
return;
