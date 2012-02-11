OBJDIR=obj
SRCDIR=src

CC = cc
CCFLAGS = -O2 -pedantic -Wall
LDFLAGS = -lrt

_OBJ = matrix.o seq.o fasta.o align.o sw.o
OBJ = $(addprefix $(OBJDIR)/,$(_OBJ))

default : sw 

$(OBJDIR) :
	mkdir $(OBJDIR)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CCFLAGS) -c -o $@ $<

sw : $(OBJDIR) $(OBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o sw $(OBJ)

.PHONY: clean

clean:
	@rm -rf $(OBJDIR)
	@rm -f sw
