## create background for HOMER by extending upstream downstream

rm(list=ls())
#fpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/'
#fpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMS/'
#outpath = '/home/bsharmi6/NA_TF_project/scMethylome/TF_first_DMR/'

## list files
#filelist = list.files(fpath)

## read each file
#for (ifile in 1:length(filelist)){
## read the TF peak/motif/DMR input file
dat = read.table('/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/mSst-2/maf/maf.txt', h=T)
#dat = out.dat[, c(1:3)]
## TF name
#TF = filelist[ifile]
dat.length = dat$end - dat$start
## upstream end
up_end = dat$start - 2000
## upstream start
up_start = up_end - dat.length
## upstream df
up_df = data.frame(chrom = dat$chrom, start = up_start, end = up_end)
## downstream start
down_start = dat$end + 2000
## downstream end
down_end = down_start + dat.length
## downstream df
down_df = data.frame(chrom = dat$chrom, start = down_start, end = down_end)
## final background df
dat_background = rbind.data.frame(up_df, down_df)
## write background and peaks to file
#write.table(dat, file = paste0(fpath, filelist[ifile], '/', filelist[ifile], '_dms.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = F)
#write.table(dat_background, file = paste0(fpath, filelist[ifile], '/', filelist[ifile], '_background.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = F)
#}
