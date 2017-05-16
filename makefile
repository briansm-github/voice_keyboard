OBJS  = vk.o 

CC=gcc 
CFLAGS= -O2
LDFLAGS=  -O -lm

all:    $(OBJS)
	$(CC) -o vk vk.o -lm  -lasound

clean: 
	rm -f *.o *~

.c.o:	$<
	$(CC) $(CFLAGS) -c $<

