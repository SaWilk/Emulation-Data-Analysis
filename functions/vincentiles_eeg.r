# vincentiles_eeg.r
#
# Computes nb_bins vincentiles for the rt vector 'dat' and corresponding eegdata. Sorts and
# averages the eegdata according to the bins. eegdata shape must (!) be timepoint x trials
# 
# Originally written by:  Trisha van Zandt
# input: 
# dat: rt (vector)
# eegdata: EEG data array (frames,trials)
# nb_bins: number of bins to compute
# adapted by Sven Hoffmann for EEG-vincentiles and translated to GNU R 15/07/2013

vincentiles_eeg <- function(dat,eegdata,nb_bins){
  # sort data and get length and dimensions
  sorted <- sort(dat,index=TRUE)
  dat <-  sorted$x
  index <-  sorted$ix
  ntrialseeg <- ncol(eegdata)
  nframeseeg <- nrow(eegdata)
  # sort eegdata according to rts
  tmpeeg <- matrix(0,nframeseeg,ntrialseeg)
  for (ind in c(1: length(index))){
    tmpeeg[,ind] <- eegdata[,index[ind]]  
  }
  
  # create temporary arrays with nb_bins replications of each value (vincentiles)
  #rts
  copyx <- rep(dat,each = nb_bins)
  vinx<-array(0,c(nb_bins,ntrialseeg))
  counter <- 1
  for (ind in c(1:nb_bins)){
    vinx[ind,] <- copyx[c(counter:(counter +ntrialseeg-1))]
    counter <-  counter + ntrialseeg
  }
  
  #eeg (eeg vincentiles according to corresponding rts)
  copyxeeg <- tmpeeg[, rep(1:ncol(tmpeeg), each = nb_bins)]
  vinxeeg<-array(0,c(nframeseeg,ntrialseeg,nb_bins))
  counter <- 1
  for (ind in c(1:nb_bins)){
    vinxeeg[,,ind] <- copyxeeg[,c(counter:(counter +ntrialseeg-1))]
    counter <-  counter + ntrialseeg
  }
  
  # obtain mean values for each (RT EEG) bin
  
  return(list(vince_eeg=apply(vinxeeg, c(1,3), mean), vince_rts=rowMeans(vinx)))
  
}