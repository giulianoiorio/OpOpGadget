CFLAGS= -O2 -Wall -fPIC -fopenmp -shared
CC= /usr/local/bin/gcc-5

EXEC_G= GenerateModel.so
EXEC_M= GeneralModel.so

OBJS_G=MT_random.o GenerateModel.o
OBJS_M= GeneralModel.o



INCL= Makefile


all: $(EXEC_M) $(EXEC_G)

$(EXEC_M): $(OBJS_M)
	$(CC) $(CFLAGS) $(OBJS_M) -o $(EXEC_M)

$(EXEC_G): $(OBJS_G)
	$(CC) $(CFLAGS) $(OBJS_G) -o $(EXEC_G)



$(OBJS_M): $(INCL)
$(OBJS_G): $(INCL)

