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
	survfitCI.Rnw \
	survfitms.Rnw \
	tail
#	coxdetail.nw

SFUN =  agsurv.R \
	coxph.R \
	predict.coxph.R \
	pyears.R \
	residuals.survreg.S\
	survConcordance.R \
	survConcordance.fit.R \
	survexp.R \
	survfitCI.R \
	survfit.coxph.R \
	survfitcoxph.fit.R \
	survfitms.R

CFUN =  agsurv4.c agsurv5.c concordance1.c coxcount1.c \
        coxexact.c \
#	msurv.c
#	coxdetail2.c

RDIR = ../R
RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)
DOCDIR= ../inst/doc

all:  doc fun

doc: $(DOCDIR)/sourcecode.pdf

$(DOCDIR)/sourcecode.pdf: all.pdf
	cp all.pdf $(DOCDIR)/sourcecode.pdf

all.pdf: all.tex
	pdflatex all   #twice to get the table of contents right
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

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
	echo "source('noweb.R'); print(notangle(nwread('all.nw'), target='$(*F)'), file='zz')" | $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
	cat zz >> $@
	rm zz
#	notangle -R$(*F) all.nw >> $@

%.S: 
	echo "# Automatically generated from all.nw using noweb" > $@
	echo "source('noweb.R'); print(notangle(nwread('all.nw'), target='$(*F)'), file='zz')" | $(R_HOME)/bin$(R_ARCH_BIN)/R   --vanilla --slave
	cat zz >> $@
	rm zz
#	notangle -R$(*F) all.nw >> $@

%.c: 
	echo "/* Automatically generated from all.nw using noweb */" > $@
	echo "source('noweb.R'); print(notangle(nwread('all.nw'), target='$(*F)'), file='zz')" |  $(R_HOME)/bin$(R_ARCH_BIN)/R  --vanilla --slave
	cat zz >> $@
	rm zz
#	notangle -R$(*F) all.nw >> $@

clean:
	-rm all.nw all.log all.aux all.toc all.tex 
	-rm $(RFUN)
	-rm $(CFUN2)
