#!/usr/bin/env python
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
