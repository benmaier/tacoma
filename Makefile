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
	/usr/local/bin/python2.7 setup.py install --user

groot:
	git fetch
	git pull
	make clean
	make grootinstall
