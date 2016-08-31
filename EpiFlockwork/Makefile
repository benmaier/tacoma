TARGET = test

default: Utilities.cpp Utilities.h Events.h Events.cpp SIS.cpp SIS.h
	g++ -O3 -std=c++11 -stdlib=libc++ main.cpp SIR.cpp SIS.cpp Events.cpp Utilities.cpp SIRS.cpp -o $(TARGET)

test_shuffle: Utilities.cpp Utilities.h Events.h Events.cpp SIS.cpp SIS.h
	g++ -O3 -std=c++11 -stdlib=libc++ test_shuffle.cpp Events.cpp Utilities.cpp -o $(TARGET)

clean:
	-rm -f *.o
	-rm -f $(TARGET)
