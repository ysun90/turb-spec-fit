function Spec = ind2spec(index,database_index)
    
    % Convert the index
    choc = database_index.index_choc(index);
    channel = database_index.index_channel(index);
    plateau = database_index.index_plateau(index);
    
    % Load spctrum data for one shot
    Spectra = load(['Spectra@',num2str(choc),'.mat']);

    % Spectrum from channel 1 or 2
    if channel==1
      Spec = Spectra.S(:,plateau); 
    elseif channel==2
      Spec = Spectra.S_2(:,plateau);
    end

end