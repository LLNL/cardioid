anat2pio: anat2pio.cc heartIO.cc EndianSwap.cc
	c++ -g -o $@ $^

writeIBMAnatomy: writeIBMAnatomy.cc
	g++ -g -o $@ $^

testReader: testReader.cc
	g++ -g -o $@ $^

clean:
	rm -f anat2pio writeIBMAnatomy testReader
