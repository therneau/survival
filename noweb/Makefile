PARTS = main.Rnw \
	coxph.Rnw \
	exact.nw \
	coxsurv.Rnw \
	coxsurv2.Rnw \
	predict.coxph.Rnw \
	concordance.Rnw \
	survexp.Rnw \
	pyears.Rnw \
	residuals.survreg.Rnw \
	survfit.Rnw \
	survfitCI.Rnw \
	msurv.nw \
	survfitms.Rnw \
	plot.Rnw \
	tail
#	coxdetail.nw

SFUN =  agsurv.R \
	coxph.R \
	plot.survfit.R \
	predict.coxph.R \
	pyears.R \
	residuals.survreg.S\
	survConcordance.R \
	survConcordance.fit.R \
	survexp.R \
	survfit.R \
	survfitCI.R \
	survfit.coxph.R \
	survfitcoxph.fit.R \
	survfitms.R

CFUN =  agsurv4.c agsurv5.c concordance1.c coxcount1.c \
        coxexact.c \
	survfitci.c
#	coxdetail2.c

RDIR = ../R
RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)
DOCDIR= ../inst/doc

all:  noweb.sty doc fun

doc: $(DOCDIR)/sourcecode.pdf

$(DOCDIR)/sourcecode.pdf: all.pdf
	cp all.pdf $(DOCDIR)/sourcecode.pdf

all.pdf: all.tex noweb.sty
	texi2dvi --pdf all.tex

all.nw: $(PARTS)
	cat $(PARTS) > all.nw

all.tex: all.nw
	echo "library(noweb); noweave('all.nw')" | $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
#	noweave -delay all.nw > all.tex

$(SFUN): all.nw
$(CFUN): all.nw
$(CFUN2): all.nw
$(RFUN): all.nw

.PHONY: fun clean doc
fun: $(RFUN) $(CFUN2)

test: $(RFUN)
	echo $(RFUN)

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	echo "require(noweb); notangle('all.nw', target='$(*F)', out='zz')" | $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
	cat zz >> $@
	rm zz

%.S: 
	echo "# Automatically generated from all.nw using noweb" > $@
	echo "require(noweb); notangle('all.nw', target='$(*F)', out='zz')" | $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
	cat zz >> $@
	rm zz

%.c: 
	echo "/* Automatically generated from all.nw using noweb */" > $@
	echo "require(noweb); notangle('all.nw', target='$(*F)', out='zz')" | $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
	cat zz >> $@
	rm zz
clean:
	-rm all.nw all.log all.aux all.toc all.tex 
	-rm $(RFUN)
	-rm $(CFUN2)
	-rm noweb.sty

noweb.sty:
	echo 'library(noweb); data(noweb); cat(noweb.sty, sep="\n", file="noweb.sty")' |  $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
