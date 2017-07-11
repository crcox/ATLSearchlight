export MCC=$(MATLABDIR)/bin/mcc
export MEX=$(MATLABDIR)/bin/mex
export MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
TOP := $(shell pwd)
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
DATA=data
GLMNET=$(DEP)/glmnet
INCL= -N -p stats -I $(SRC) -I $(GLMNET)
.PHONEY: all clean-all clean-postbuild glmnet sdist

all: glmnet ATLSearchlight clean-postbuild

glmnet: $(GLMNET)/glmnetMex.mexa64
$(GLMNET)/glmnetMex.mexa64: $(GLMNET)/glmnetMex.F $(GLMNET)/GLMnet.f
	$(MEX) -fortran -outdir $(GLMNET) $^

ATLSearchlight: $(SRC)/ATLSearchlight.m
	$(MCC) -v $(MFLAGS) $(INCL) -a $(DATA) -o $@ $<

clean-postbuild:
	-rm *.dmr
	-rm mccExcludedFiles.log
	-rm readme.txt

sdist:
	tar czhf $(SRCTAR) src dependencies data

clean-all:
	-rm ATLSearchlight
	-rm $(DEP)
	-rm $(DATA)
	-rm $(SRC)
	-rm requiredMCRProducts.txt
	-rm build-csf-ATLSearchlight.sh.e*
	-rm build-csf-ATLSearchlight.sh.o*
