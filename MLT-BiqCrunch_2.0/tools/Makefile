TARGETS = kc2bc qp2bc

CC = gcc
CFLAGS = -std=c99

default : $(TARGETS)

% : %.c
	$(CC) $(CFLAGS) -o $@ $<

clean :
	$(RM) $(TARGETS) *.bc *output*

