PARTS = main.Rnw \
	residuals.survreg.Rnw \
	coxsurv.Rnw \
	survexp.Rnw \
	predict.coxph.Rnw 

SFUN =  residuals.survreg.S\
	survfit.coxph.R \
	agsurv.R\
	survexp.R \
	predict.coxph.S 

CFUN =  agsurv4.c agsurv5.c

all.pdf: all.tex
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw
	echo "\\end{document}" >> all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

RDIR = ../../R
RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../../src/%)

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
