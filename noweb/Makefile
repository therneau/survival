PARTS = main.Rnw \
	coxph.Rnw \
	exact.nw \
	agreg.Rnw \
	coxsurv.Rnw \
	coxsurv2.Rnw \
	coxsurv3.Rnw \
	finegray.Rnw \
	predict.coxph.Rnw \
	concordance.Rnw \
	survexp.Rnw \
	parse.Rnw \
	pyears.Rnw pyears2.Rnw \
	residuals.survfit.Rnw \
	residuals.survfit2.Rnw \
	residuals.survreg.Rnw \
	survfit.Rnw \
	survfitKM.Rnw \
	survfitCI.Rnw \
	msurv.nw \
	survfitms.Rnw \
	survexpm.Rnw\
	plot.Rnw \
	statefig.Rnw\
	tmerge.Rnw\
	yates.Rnw  yates2.Rnw\
	zph.Rnw \
	tail
#	coxdetail.nw
#	residuals.survfit2.Rnw\

SFUN =  agreg.fit.R \
	agsurv.R \
        concordance.R \
	coxph.R \
	coxsurvfit.R \
	finegray.R \
	model.matrix.coxph.R \
	parsecovar.R \
	plot.survfit.R \
	predict.coxph.R \
	pyears.R \
	print.pyears.R \
	residuals.survfit.R \
	residuals.survfitcox.R \
	residuals.survreg.R\
	statefig.R \
	survexp.R \
	survexpm.R \
	survfit.R \
	survfitCI.R \
	survfitKM.R \
	survfit.coxph.R \
	survfit.coxphms.R \
	survfitms.R\
	tmerge.R \
	yates.R \
	cox.zph.R

CFUN =  agfit4.c \
        agsurv4.c agsurv5.c \
	cdecomp.c \
        concordance3.c \
        coxexact.c \
	survfitci.c \
	survfitkm.c \
#	coxdetail2.c

RDIR = ../R
RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)
DOCDIR= ../inst/doc

all:  noweb.sty doc fun

doc: code.pdf 

code.pdf: code.tex noweb.sty
	pdflatex code.tex
	pdflatex code.tex

code.nw: $(PARTS)
	cat $(PARTS) > code.nw

code.tex: code.nw
	echo "library(noweb); noweave('code.nw')" | R  --slave

$(SFUN): code.nw
$(CFUN): code.nw
$(CFUN2): code.nw
$(RFUN): code.nw

.PHONY: fun clean doc all
fun: $(RFUN) $(CFUN2) noweb.sty

test: $(RFUN)
	echo $(RFUN)

%.R: 
	echo "# Automatically generated from the noweb directory" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='z$(@F)')" | R  --slave
	cat z$(@F) >> $@
	rm z$(@F)

%.S: 
	echo "# Automatically generated from the noweb directory" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='z$(@F)')" | R  --slave
	cat z$(@F) >> $@
	rm z$(@F)

%.c: 
	echo "/* Automatically generated from the noweb directory */" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='z$(@F)')" | R  --slave
	cat z$(@F) >> $@
	rm z$(@F)

clean:
	-rm -f code.nw code.log code.aux code.toc code.tex code.bbl code.blg code.out
	-rm -f noweb.sty

noweb.sty:
	echo 'library(noweb); data(noweb); cat(noweb.sty, sep="\n", file="noweb.sty")' |  R  --slave


