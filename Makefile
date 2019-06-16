PROJ = quadrature
OBJ = quadrature.o main.o

CFLAGS = -Wall -g
LDFLAGS = -lm

all: $(PROJ)

$(PROJ): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ) $(PROJ)
