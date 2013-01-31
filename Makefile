all: rseqc samtools

rseqc:
	cd rseqc; PYTHONPATH=../resources/usr/local/lib/python2.7/site-packages:$PYTHONPATH python setup.py install --prefix=../resources/usr/local

samtools:
	make -C samtools samtools
	cp -a samtools/samtools resources/usr/local/bin/samtools

clean:
	rm -rf resources/usr/local/bin/* resources/usr/local/lib/python2.7

.PHONY: all rseqc samtools clean
