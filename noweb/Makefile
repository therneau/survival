PARTS = main.Rnw \
	coxsurv.Rnw \
	coxsurv2.Rnw \
	predict.coxph.Rnw \
	survexp.Rnw \
	residuals.survreg.Rnw 

SFUN =  residuals.survreg.S\
	survfit.coxph.R \
	survfit.coxph.fit.R\
	agsurv.R\
	survexp.R \
	predict.coxph.S 

CFUN =  agsurv4.c agsurv5.c
RDIR = ../R
RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)
DOCDIR= ../inst/doc

all:  $(DOCDIR)/sourcecode.pdf fun

$(DOCDIR)/sourcecode.pdf: all.pdf
	cp all.pdf $(DOCDIR)

all.pdf: all.tex
	latex all   #twice to get the table of contents right
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw
	echo "\\end{document}" >> all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

$(SFUN): all.nw
$(CFUN): all.nw
$(CFUN2): all.nw
$(RFUN): all.nw

.PHONY: fun clean
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
