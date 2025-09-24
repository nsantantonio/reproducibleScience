#PARENTPATH=/home/nsant/Dropbox
#cd ${PARENTPATH}/reproducibleScience
Rscript cleanDataExample.R
Rscript dataAnalysis.R

cd manuscript
pdflatex papertemplate.tex 