# Makefile for building POM with netCDF output
F77 = pgf90
FFLAGS= -O2

SRC = \
main.f advave.f advt.f advq.f advu.f advv.f archive.f baropg.f bcdata2.f \
bcond2.f dens.f display.f drifter3.f extrnl.f first.f maxmin.f \
print.f prints.f profq.f proft.f profu.f profv.f prtxy.f setdom.f \
sinter.f slicexz.f sliceyz.f smag.f tands.f vertvl.f wreal.f \
zeroes.f cutpasses.f ftest0.f ftest1.f ftest2.f ftest3.f ftest4.f\
bubsort.f genran2.f rdcauseways.f \
roulette.f  mutate2.f

SRC1 = \
main1.f advave.f advt.f advq.f advu.f advv.f archive.f baropg.f bcdata2.f \
bcond2.f dens.f display.f drifter3.f extrnl.f first.f maxmin.f \
print.f prints.f profq.f proft.f profu.f profv.f prtxy.f setdom.f \
sinter.f slicexz.f sliceyz.f smag.f tands.f vertvl.f wreal.f \
zeroes.f cutpasses.f ftest0.f ftest1.f ftest2.f ftest3.f ftest4.f\
bubsort.f genran2.f rdcauseways.f \
roulette.f  mutate2.f

SRC2 = genmain.f rdcauseways.f mkrandgrid.f mksid.f genran.f roulette.f bubsort.f 

SRC3 = propagate.f rdfitness.f roulette.f genran.f bubsort.f \
mutate2.f genran2.f

OBJS = \
advave.o advt.o advq.o advu.o advv.o archive.o baropg.o bcdata2.o \
bcond2.o dens.o display.o drifter3.o main.o extrnl.o first.o maxmin.o \
print.o prints.o profq.o proft.o profu.o profv.o prtxy.o setdom.o \
sinter.o slicexz.o sliceyz.o smag.o tands.o vertvl.o wreal.o \
zeroes.o cutpasses.o

INC = \
comdeck Makefile

DATA = \
run_data4 rdtime.dat seedhour.dat tands.dat model_grid 

otbpt:  $(SRC) genetic comdeck
	$(F77) $(FFLAGS) $(SRC) -o otbpt

genmain:$(SRC2) genetic comdeck
	$(F77) $(FFLAGS) $(SRC2) -o genmain

propagate:$(SRC3) genetic 
	$(F77) $(FFLAGS) $(SRC3) -o propagate

run_one:$(SRC1) comdeck
	$(F77) $(FFLAGS) $(SRC1) -o run_one
clean:
	rm -f ecomtest $(OBJS)

tar:
	tar hcvf - $(SRC) $(INC) $(DATA) | compress > ecomtest.tar

zip:
	zip -l ecomtest.zip $(SRC) $(INC) $(DATA) 

rtest:  rtest.f
	$(F77) $(FFLAGS) rtest.f -o rtest

temp:	temp.f
	$(F77) $(FFLAGS) temp.f -o temp
