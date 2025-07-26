SHELL := /bin/bash
# Set dynamic library flag
UNAMES=$(shell uname -s)
ifeq ($(UNAMES),Linux)
	DLIBFLAG=-shared
    PYTHONINC := $(shell python3-config --includes)
    PYTHONLIB := $(shell python3-config --ldflags)
endif
ifeq ($(UNAMES),Darwin)
	DLIBFLAG=-dynamiclib
	PYTHONINC := $(shell python3-config --includes)
    PYTHONLIB := $(shell python3-config --ldflags)
endif


XSAT_ := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
OUT_=$(XSAT_)/out
R_SQUARE_=$(OUT_)/R_square
R_ULP_=$(OUT_)/R_ulp
R_VERIFY_=$(OUT_)/R_verify

XSAT_GEN=$(XSAT_)/xsat_gen.py

ifdef IN
   $(shell echo $(IN) > XSAT_IN.txt)
endif


IN:= $(shell cat XSAT_IN.txt)

PYTHON_H:=

define XSAT_echo
	@echo "[XSat] $1 "
endef


all:  compile

gen:  build/foo.c xsat_gen.py
build/foo.c: $(IN)  XSAT_IN.txt
	@echo "[XSAT] .smt2 -> .c"
	@mkdir -p build
	python xsat_gen.py $<  > $@

compile_square: build/R_square/foo.so
build/R_square/foo.so: build/foo.c include/R_square/xsat.h  $(IN)
	@echo [XSAT]Compiling the representing function as $@
	@mkdir -p build/R_square
	@clang -O3 -fPIC $< $(DLIBFLAG) -o $@  $(PYTHONINC) -I include/R_ulp  $(PYTHONLIB) -fbracket-depth=3000

compile_verify: build/R_verify/foo.so
build/R_verify/foo.so: build/foo.c include/R_verify/xsat.h  $(IN)
	@echo [XSAT]Compiling the representing function as $@
	@mkdir -p build/R_verify
	@clang -O3 -fPIC $< $(DLIBFLAG) -o $@  $(PYTHONINC) -I include/R_ulp  $(PYTHONLIB) -fbracket-depth=3000

compile_ulp: build/R_ulp/foo.so
build/R_ulp/foo.so: build/foo.c include/R_ulp/xsat.h  $(IN)
	@echo [XSAT]Compiling the representing function as $@
	@mkdir -p build/R_ulp
	@clang -O3 -fPIC $< $(DLIBFLAG) -o $@  $(PYTHONINC) -I include/R_ulp  $(PYTHONLIB) -fbracket-depth=3000

compile:  compile_ulp compile_verify

solve: compile
	@echo [XSAT] Executing the solver.
	@python xsat.py

test: test_benchmarks.py
	python $<

helloworld: Benchmarks/div3.c.50.smt2
	make IN=$>
	python xsat.py

clean:
	$(XSAT_echo) Cleaning build/ and Results/
	@rm -vf build/foo.c build/foo.symbolTable
	@rm -vfr build/R_square build/R_ulp build/R_verify
	@rm -vf Results/*


.PHONY: copy gen clean compile compile_square test


