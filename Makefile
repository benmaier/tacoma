TARGET = test
EPIFOL = ./cFlockwork

default: $(EPIFOL)/Utilities.cpp $(EPIFOL)/Utilities.h $(EPIFOL)/Events.h $(EPIFOL)/Events.cpp $(EPIFOL)/SIS.cpp $(EPIFOL)/SIS.h $(EPIFOL)/SIR.h $(EPIFOL)/SIR.cpp $(EPIFOL)/SIRS.cpp $(EPIFOL)/SIRS.h
	g++ -O3 -std=c++11 -stdlib=libc++ -I$(EPIFOL) $(EPIFOL)/main.cpp $(EPIFOL)/SIR.cpp $(EPIFOL)/SIS.cpp $(EPIFOL)/Events.cpp $(EPIFOL)/Utilities.cpp $(EPIFOL)/SIRS.cpp $(EPIFOL)/EqFlockwork.cpp $(EPIFOL)/SI.cpp -o $(TARGET)

test_shuffle: $(EPIFOL)/Utilities.cpp $(EPIFOL)/Utilities.h $(EPIFOL)/Events.h $(EPIFOL)/Events.cpp $(EPIFOL)/SIS.cpp $(EPIFOL)/SIS.h
	g++ -O3 -std=c++11 -stdlib=libc++ -I$(EPIFOL) $(EPIFOL)/test_shuffle.cpp $(EPIFOL)/Events.cpp $(EPIFOL)/Utilities.cpp -o $(TARGET)

clean:
	-rm -f *.o
	-rm -f $(TARGET)

clean_all:
	make clean
	make pyclean
	make matclean

pyclean:
	-rm -f *.so
	-rm -rf *.egg-info*
	-rm -rf ./tmp/
	-rm -rf ./build/

matclean:
	-rm -rf ./matlabbuild/
