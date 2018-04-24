default: 
	make python

clean:
	-rm -f *.o
	make pyclean

clean_all:
	make clean
	make pyclean

pyclean:
	-rm -f *.so
	-rm -rf *.egg-info*
	-rm -rf ./tmp/
	-rm -rf ./build/

python:
	pip install -e ../tacoma --no-binary :all:

grootinstall:
	/usr/local/bin/pip2.7 install --user ../tacoma

groot:
	git fetch
	git pull
	make clean
	make grootinstall
