.PHONY:all compile test clean

all: compile test

compile:
	(cd src; make)
test: compile
	(cd test; make)

clean:
	(cd src; make clean)
	(cd test; make clean)

