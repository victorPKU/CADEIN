SHELL           = /bin/sh

CC              = gcc

CC_FLAGS        = -ansi -O6 -fomit-frame-pointer -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing

CC_LINKERS      = -lm

STRIP           = strip

SECURITY	= chmod 555


CC_FLAGS_FULL	=  $(CC_FLAGS)


.SUFFIXES:	.c .o

.c.o:
		$(CC) $(CC_FLAGS_FULL) -c $<  


PROGRAMS = cadein buildmod

all:		$(PROGRAMS)
 

cadein:		cadein.o datatype.o rotamer.o record.o mem.o coordination.o combine.o protein.o mutation.o oxypair.o oxylib.o backrub.o geometry.o superpose.o mem.h rotamer.h record.h combine.h oxypair.h oxylib.h mutation.h backrub.h protein.h
		$(CC) $(CC_FLAGS) -o $@ cadein.o datatype.o rotamer.o record.o mem.o coordination.o combine.o protein.o mutation.o oxypair.o oxylib.o backrub.o geometry.o superpose.o $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@


buildmod:		buildmod.o datatype.o mem.o record.o mutation.o geometry.o rotamer.o coordination.o combine.o protein.o backrub.o superpose.o model.h record.h mutation.h geometry.h rotamer.h mem.h backrub.h protein.h
		$(CC) $(CC_FLAGS) -o $@ buildmod.o datatype.o mem.o record.o mutation.o geometry.o rotamer.o coordination.o combine.o protein.o backrub.o superpose.o $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@


clean:
		rm -f *.o core $(PROGRAMS)


cadein.o:			mem.h rotamer.h record.h combine.h oxypair.h oxylib.h mutation.h backrub.h protein.h
datatype.o:              datatype.h
mem.o:			mem.h
buildmod.o:              model.h record.h mutation.h geometry.h rotamer.h mem.h backrub.h protein.h
record.o:			record.h mem.h mutation.h
mutation.o:			mutation.h mem.h geometry.h rotamer.h superpose.h
geometry.o:			geometry.h 
rotamer.o:			rotamer.h
superpose.o:             superpose.h geometry.h mem.h
backrub.o:               backrub.h geometry.h mem.h
combine.o:               combine.h coordination.h mem.h 
coordination.o:          coordination.h geometry.h 
oxypair.o:               oxypair.h mem.h geometry.h
oxylib.o:                oxylib.h mem.h geometry.h
protein.o:               protein.h mem.h geometry.h 


