OBJDIR=obj
SRCDIR=src

CC = cc
CFLAGS = -O2 -ansi -pedantic -Wall
LDFLAGS =

_OBJ = matrix.o seq.o fasta.o align.o sw.o
OBJ = $(addprefix $(OBJDIR)/,$(_OBJ))

default : sw 

$(OBJDIR) :
	mkdir $(OBJDIR)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

sw : $(OBJDIR) $(OBJ)
	$(CC) -o sw $(OBJ) $(LDFLAGS)

.PHONY: clean

clean:
	@rm -rf $(OBJDIR)
	@rm -f sw
