FLAGS=-O2 \
	 -fopenmp -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -ftree-vectorize -fipa-cp-clone -w

CC=gcc

RM=rm -f

EXEC=dna

all: $(EXEC)

$(EXEC):
	$(CC) $(FLAGS) dna.c -c -o dna.o
	$(CC) $(FLAGS) dna.o -o $(EXEC)

run:
	./$(EXEC)

clean:
	$(RM) dna.o $(EXEC)