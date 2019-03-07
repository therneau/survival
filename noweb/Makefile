PARTS = main.Rnw \
	coxph.Rnw \
	exact.nw \
	agreg.Rnw \
	coxsurv.Rnw \
	coxsurv2.Rnw \
	finegray.Rnw \
	predict.coxph.Rnw \
	concordance.Rnw \
        survconcordance.Rnw \
	survexp.Rnw \
	parse.Rnw \
	pyears.Rnw pyears2.Rnw \
	residuals.survreg.Rnw \
	survfit.Rnw \
	survfitCI.Rnw \
	survfitKM.Rnw \
	msurv.nw \
	survfitms.Rnw \
	plot.Rnw \
	statefig.Rnw\
	tmerge.Rnw\
	yates.Rnw  yates2.Rnw\
	tail
#	coxdetail.nw

SFUN =  agreg.fit.R \
	agsurv.R \
        concordance.R \
	coxph.R \
	finegray.R \
	model.matrix.coxph.R \
	parsecovar.R \
	plot.survfit.R \
	predict.coxph.R \
	pyears.R \
	print.pyears.R \
	residuals.survreg.R\
	statefig.R \
	survConcordance.R \
	survConcordance.fit.R \
	survexp.R \
	survfit.R \
	survfitCI.R \
	survfitKM.R \
	survfit.coxph.R \
	survfitcoxph.fit.R \
	survfitms.R\
	tmerge.R \
	yates.R

CFUN =  agfit4.c \
        agsurv4.c agsurv5.c \
        concordance1.c concordance3.c coxcount1.c \
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
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='zz')" | R  --slave
	cat zz >> $@
	rm zz

%.S: 
	echo "# Automatically generated from the noweb directory" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='zz')" | R  --slave
	cat zz >> $@
	rm zz

%.c: 
	echo "/* Automatically generated from the noweb directory */" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='zz')" | R  --slave
	cat zz >> $@
	rm zz

clean:
	-rm code.nw code.log code.aux code.toc code.tex code.bbl code.blg code.out
	-rm noweb.sty

noweb.sty:
	echo 'library(noweb); data(noweb); cat(noweb.sty, sep="\n", file="noweb.sty")' |  R  --slave
