FLAGS=-O3 -Wunused-result

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