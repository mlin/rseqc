# This file is a part of rseqc (DNAnexus platform app).
# Copyright (C) 2013 DNAnexus, Inc.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

all: rseqc samtools

rseqc:
	cd rseqc; PYTHONPATH=../resources/usr/local/lib/python2.7/site-packages:$PYTHONPATH python setup.py install --prefix=../resources/usr/local

samtools:
	make -C samtools samtools
	cp -a samtools/samtools resources/usr/local/bin/samtools

clean:
	$(MAKE) -C samtools clean
	rm -rf resources/usr/local/bin/* resources/usr/local/lib/python2.7/site-packages/*

.PHONY: all rseqc samtools clean
