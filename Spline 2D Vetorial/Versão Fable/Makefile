# Compilação dos dois programas (gfortran + MKL). Uso:
#   source /opt/intel/oneapi/setvars.sh && make          (MKL, padrão)
#   make LAPACK=sistema                                  (LAPACK/BLAS do sistema)
FC      = gfortran
FFLAGS  = -O3 -fopenmp -m64
ifeq ($(LAPACK),sistema)
  LIBS  = -llapack -lblas
else
  LIBS  = -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif
SRC     = src
MAIN    = $(SRC)/2D_Spline_Vetorial_Diferentes_zeta_omp.f90
PLOT    = $(SRC)/Plot_funcoes_zeta.f90
INCS    = $(wildcard $(SRC)/kernel.f90 $(SRC)/gamma_min.f90)

all: bin/spline_omp bin/plot

bin/spline_omp: $(MAIN) $(INCS)
	@mkdir -p bin
	cd $(SRC) && $(FC) $(FFLAGS) -J ../bin $(notdir $(MAIN)) -o ../bin/spline_omp $(LIBS)
	@echo "$(FC) $(FFLAGS) $(LIBS)" > bin/flags.txt

bin/plot: $(PLOT) $(INCS)
	@mkdir -p bin
	cd $(SRC) && $(FC) -O2 -J ../bin $(notdir $(PLOT)) -o ../bin/plot

clean:
	rm -f bin/spline_omp bin/plot bin/*.mod bin/flags.txt

.PHONY: all clean
