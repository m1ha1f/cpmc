#flex -P__cfl_ -s -o$@ $<

CXX = g++ -g -O3 -Wall
IFLAGS = -Iutil -Igroup -Ipartition 
LFLAGS = -lm -ljpeg
F90FLAGS = -lgfortran -llapack -lblas 

#######################################################################################333

UTILOBJS = util/configFileLexer.o util/configure.o util/exception.o util/random.o util/regex.o \
						util/string.o util/timer.o util/util.o util/kmeans.o util/image.o util/message.o util/segmentation.o

GROUPOBJS = group/filterbank.o group/nonmax.o group/histogram.o group/texture.o \
							group/color.o group/pb.o group/ic.o group/region.o group/affinity.o

PARTOBJS = partition/smatrix.o partition/dmatrix.o

#######################################################################################333

.SUFFIXES : .o .cc 

.cc.o:
	${CXX} -c ${IFLAGS} -o $@ -x c++ $<

#######################################################################################333

all: pixelfeatures unitex segment textons oe pr_calc segdump

#######################################################################################333

segment: util/libutil.a group/libgroup.a partition/libpartition.a trlan/libtrlan.a apps/segment.o
	${CXX} -v -o $@ apps/segment.o group/libgroup.a util/libutil.a partition/libpartition.a \
	trlan/libtrlan.a ${IFLAGS} ${LFLAGS} ${F90FLAGS}

segdump: util/libutil.a group/libgroup.a partition/libpartition.a trlan/libtrlan.a apps/segdump.o
	${CXX} -v -o $@ apps/segdump.o group/libgroup.a util/libutil.a partition/libpartition.a \
	trlan/libtrlan.a ${IFLAGS} ${LFLAGS} ${F90FLAGS} 

pixelfeatures: util/libutil.a group/libgroup.a apps/pixelfeatures.o
	${CXX} -v -o $@ apps/pixelfeatures.o group/libgroup.a util/libutil.a ${IFLAGS} ${LFLAGS}

oe: util/libutil.a group/libgroup.a apps/oe.o
	${CXX} -v -o $@ apps/oe.o group/libgroup.a util/libutil.a ${IFLAGS} ${LFLAGS}

unitex: util/libutil.a group/libgroup.a apps/unitex.o
	${CXX} -v -o $@ apps/unitex.o group/libgroup.a util/libutil.a ${IFLAGS} ${LFLAGS}

textons: util/libutil.a group/libgroup.a apps/textons.o
	${CXX} -v -o $@ apps/textons.o group/libgroup.a util/libutil.a ${IFLAGS} ${LFLAGS}

pr_calc: util/libutil.a group/libgroup.a partition/libpartition.a bench/match.o bench/pr_calc.o csa
	${CXX} -v -o $@ bench/pr_calc.o bench/match.o group/libgroup.a util/libutil.a

#######################################################################################333

util/libutil.a: ${UTILOBJS} 
	rm -f util/libutil.a && ar cr util/libutil.a $(UTILOBJS) && ranlib util/libutil.a

group/libgroup.a: ${GROUPOBJS} 
	rm -f group/libgroup.a && ar cr group/libgroup.a $(GROUPOBJS) && ranlib group/libgroup.a

partition/libpartition.a: ${PARTOBJS} 
	rm -f partition/libpartition.a && ar cr partition/libpartition.a $(PARTOBJS) && ranlib partition/libpartition.a

trlan/libtrlan.a:
	+(cd trlan && make -f Makefile.gcc)

csa:
	+(cd bench/csa-1.2/csa/prec_costs && make install)

#######################################################################################333

clean :
	rm -f group/*.a util/*.a partition/*.a trlan/*.a 
	rm -f apps/*.o bench/*.o group/*.o util/*.o partition/*.o 
	rm -f apps/*~ group/*~ util/*~ partition/*~

cleaner : clean
	rm -f segment segdump pixelfeatures unitex textons pr_calc oe 

cleanest : cleaner
	(cd trlan && make -f Makefile.gcc clean)
	(cd bench/csa-1.2/csa/prec_costs && make clean)
	rm -f csa
