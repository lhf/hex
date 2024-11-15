R= 4
I= -
O= out

all:
	python hex.py $R $I > $O.eps
	-pstopdf $O.eps
	grep '^%=CSV ' $O.eps | cut -c7- > $O.csv

clean:
	-rm -f $O.*

.PHONY:	all clean
