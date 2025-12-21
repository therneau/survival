# Over the next few years the plan is to migrate all of the material out of
#  noweb, with documentation in the methods vignette. I like noweb (T Therneau)
#  but its set of users is near zero, and the package needs to move to a more
#  open model of maintainance and developement.
#
PARTS = main.Rnw \
	coxph.Rnw \
	exact.nw \
	agreg.Rnw \
	coxsurv.Rnw \
	coxsurv3.Rnw \
	finegray.Rnw \
	predict.coxph.Rnw \
	survexp.Rnw \
	parse.Rnw \
	pyears.Rnw pyears2.Rnw \
	residuals.survfit.Rnw \
	residuals.survfit2.Rnw \
	residuals.survreg.Rnw \
	survfit.Rnw \
	msurv.nw \
	statefig.Rnw\
	yates.Rnw  yates2.Rnw\
	zph.Rnw \
	tail
#	concordance.Rnw \
#	coxsurv2.Rnw \
#	coxdetail.nw
#	plot.Rnw \
#	residuals.survfit2.Rnw\
#	tmerge.Rnw\

SFUN =  agreg.fit.R \
	coxph.R \
	finegray.R \
	parsecovar.R \
	predict.coxph.R \
	pyears.R \
	print.pyears.R \
	residuals.survreg.R\
	statefig.R \
	survexp.R \
	survfit.coxphms.R \
	yates.R \
	cox.zph.R
#       concordance.R \
#	model.matrix.coxph.R \
#	plot.survfit.R \
#	residuals.survfitcox.R \
#	survfit.R \
#	tmerge.R \

# CFUN =  concordance3.c \
CFUN =  coxexact.c

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


