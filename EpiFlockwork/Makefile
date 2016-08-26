TARGET = SIStest

default: Utilities.cpp Utilities.h Events.h Events.cpp SIS.cpp SIS.h
	g++ -lm -O3 -std=c++11 -stdlib=libc++ main.cpp SIS.cpp Events.cpp Utilities.cpp -o $(TARGET)

clean:
	-rm -f *.o
	-rm -f $(TARGET)
