#!/usr/bin/env python
#
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

import os, sys, unittest, json, subprocess

import dxpy, dxpy.app_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")
preview = 50

class TestMyApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        bundled_resources = dxpy.app_builder.upload_resources(src_dir)
        program_id = dxpy.app_builder.upload_applet(src_dir, bundled_resources, overwrite=True)
        cls.program = dxpy.DXApplet(program_id)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_unpaired(self):
        bed_file = dxpy.find_one_data_object(name="hg19_GRCh37_Feb2009_RefSeq.bed")['id']
        mappings = dxpy.find_one_data_object(name="unpaired_RNA-Seq_mappings", typename="LetterMappings")['id']
        if bed_file == None:
            print "Cannot find hg19_GRCh37_Feb2009_RefSeq.bed.  Please upload it"
            return False
        if mappings == None:
            print "Cannot find unpaired_RNA-Seq_mappings.  Please upload it"
            return False

        input = { 'rna_seq_mappings': dxpy.dxlink(mappings), 
                  'bed_file': dxpy.dxlink(bed_file) }

        print "Running program with", input
        job = self.program.run(input)
        print "launched test_unpaired ", job.get_id()

    def test_small_paired_with_contam(self):
        bed_file = dxpy.find_one_data_object(name="hg19_GRCh37_Feb2009_RefSeq.bed")['id']
        mappings = dxpy.find_one_data_object(name="SRR018256_end100000_mappings", typename="LetterMappings")['id']
        contam_contig = dxpy.find_one_data_object(name="human rRNA", typename="ContigSet")['id']
        reads = dxpy.find_one_data_object(name="SRR018256_1_end100000 reads", typename="LetterReads")['id']
        if bed_file == None:
            print "Cannot find hg19_GRCh37_Feb2009_RefSeq.bed. Please upload it."
            return False
        if mappings == None:
            print "Cannot find unpaired_RNA-Seq_mappings. Please upload it."
            return False
        if contam_contig == None:
            print "can't find contigset"
            return False
        if reads == None:
            print "can't find reads"
            return False

        input = { 'rna_seq_mappings': dxpy.dxlink(mappings), 
                  'bed_file': dxpy.dxlink(bed_file),
                  'contaminants': [dxpy.dxlink(contam_contig)],
                  'original_reads': [dxpy.dxlink(reads)] }

        print "Running program with", input
        job = self.program.run(input)
        print "launched test_unpaired ", job.get_id()

    def test_paired_with_contam(self):
        bed_file = dxpy.find_one_data_object(name="hg19_GRCh37_Feb2009_RefSeq.bed")['id']
        mappings = dxpy.find_one_data_object(name="SRR018256_paired_RNA_Mappings", typename="LetterMappings")['id']
        contam_contig = dxpy.find_one_data_object(name="human rRNA", typename="ContigSet")['id']
        reads = dxpy.find_one_data_object(name="SRR018256_reads", typename="LetterReads")['id']
        if bed_file == None:
            print "Cannot find hg19_GRCh37_Feb2009_RefSeq.bed.  Please upload it"
            return False
        if mappings == None:
            print "Cannot find Mappings.  Please upload them"
            return False
        if contam_contig == None:
            print "Cannot find human rRNA.  Please upload it"
            return False
        if reads == None:
            print "Cannot find SRR018256_reads.  Please upload it"
            return False

        input = { 'rna_seq_mappings': dxpy.dxlink(mappings), 
                  'bed_file': dxpy.dxlink(bed_file),
                  'contaminants': [dxpy.dxlink(contam_contig)],
                  'original_reads': [dxpy.dxlink(reads)] }

        print "Running program with", input
        job = self.program.run(input)
        print "launched test_paired_with_contam ", job.get_id()
        
if __name__ == '__main__':
    unittest.main()
