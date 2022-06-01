function OUTEEG = dfgremICs2021(INEEG,save)
    TMPEEG = INEEG;
    
    TMPEEG = iclabel(TMPEEG);
    artifacts  = find(TMPEEG.etc.ic_classification.ICLabel.classifications(:,1) > 0.6);
    % remove ocular artifacts
    TMPEEG = pop_subcomp( TMPEEG, artifacts,0);

    % save data if desired.
    if save == true 
       TMPEEG=pop_saveset(TMPEEG,'filename',[TMPEEG.filename(1:end-4) '-clean'],...
        'filepath',TMPEEG.filepath);
    end
    OUTEEG = TMPEEG;
    return
    