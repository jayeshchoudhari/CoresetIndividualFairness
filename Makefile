all: mainIF_MV mainIF_Coreset

mainIF_MV: mainIF_MV.o dataIO.o util.o
	g++ -O3 mainIF_MV.o dataIO.o util.o -o mainIF_MV

mainIF_Coreset: mainIF_Coreset.o dataIO.o util.o
	g++ -O3 mainIF_Coreset.o dataIO.o util.o -o mainIF_Coreset


mainIF_MV.o: mainIF_MV.cpp ./include/dataIO.cpp  ./include/dataIO.h ./include/util.cpp  ./include/util.h ./include/namespace.h
	g++ -O3 -c mainIF_MV.cpp

mainIF_Coreset.o: mainIF_Coreset.cpp ./include/dataIO.cpp  ./include/dataIO.h ./include/util.cpp  ./include/util.h ./include/namespace.h
	g++ -O3 -c mainIF_Coreset.cpp

dataIO.o: ./include/dataIO.cpp  ./include/dataIO.h ./include/util.cpp  ./include/util.h ./include/namespace.h
	g++ -O3 -c ./include/dataIO.cpp

util.o: ./include/util.cpp  ./include/util.h  ./include/namespace.h
	g++ -O3 -c ./include/util.cpp

clean:
	rm *.o mainIF_MV mainIF_Coreset