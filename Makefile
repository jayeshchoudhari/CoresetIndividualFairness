mainIF_MV: mainIF_MV.o dataIO.o util.o
	g++ -O3 mainIF_MV.o dataIO.o util.o -o mainIF_MV

mainIF_MV.o: ./include/dataIO.cpp  ./include/dataIO.h ./include/util.cpp  ./include/util.h ./include/namespace.h
	g++ -O3 -c mainIF_MV.cpp

dataIO.o: ./include/dataIO.cpp  ./include/dataIO.h ./include/util.cpp  ./include/util.h ./include/namespace.h
	g++ -O3 -c ./include/dataIO.cpp

util.o: ./include/util.cpp  ./include/util.h  ./include/namespace.h
	g++ -O3 -c ./include/util.cpp

clean:
	rm *.o mainIF_MV