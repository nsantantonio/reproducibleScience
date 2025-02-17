#PARENTPATH=/home/nsant/Dropbox
#cd ${PARENTPATH}/reproducibleScience
Rscript cleanData.R
Rscript dataAnalysis.R

cd manuscript
pdflatex papertemplate.tex 