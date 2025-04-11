CC = cc
RCMD  = R
RCMDARGS = CMD SHLIB
RM = rm
SRCs =	RadonTrans.c Field2D.c \
	gaussrand.c poissonrand.c mtrand.c

all: RSHLIB

clean:
	$(RM) -f *.o *.so a.out


RSHLIB:
	$(RCMD) $(RCMDARGS) $(SRCs)

TEST:
	$(CC) -D_TEST_RADONTRANS_ $(SRCs)
