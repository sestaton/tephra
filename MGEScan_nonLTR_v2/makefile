CC=     gcc
CFLAG= -O3 -Wall
SRCS=   translate.c
OBJS=   translate.o
EXES=   translate

$(OBJS):	$(SRCS)   
	$(CC) $(CFLAG) -c $(SRCS)
all:	$(OBJS)
	$(CC) $(CFLAG) -o $(EXES) $(OBJS) 
clean:	
	rm -rf *.o translate *~
