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

        input = { 'RNA-Seq Mappings': dxpy.dxlink(mappings), 
                  'BED file': dxpy.dxlink(bed_file) }

        print "Running program with", input
        job = self.program.run(input)
        print "Waiting for test_unpaired ", job.get_id()
        job.wait_on_done()
        print json.dumps(job.describe()["output"])


    def test_unpaired_with_contam(self):
        bed_file = dxpy.find_one_data_object(name="hg19_GRCh37_Feb2009_RefSeq.bed")['id']
        mappings = dxpy.find_one_data_object(name="unpaired_RNA-Seq_mappings", typename="LetterMappings")['id']
        contam_contig = dxpy.find_one_data_object(name="human rRNA", typename="ContigSet")['id']
        reads = dxpy.find_one_data_object(name="SRR399297 reads", typename="LetterReads")['id']
        if bed_file == None:
            print "Cannot find hg19_GRCh37_Feb2009_RefSeq.bed.  Please upload it"
            return False
        if mappings == None:
            print "Cannot find unpaired_RNA-Seq_mappings.  Please upload it"
            return False
        if contam_contig == None:
            print "Cannot find human rRNA.  Please upload it"
            return False
        if reads == None:
            print "Cannot find SRR399297 reads.  Please upload it"
            return False


        input = { 'RNA-Seq Mappings': dxpy.dxlink(mappings), 
                  'BED file': dxpy.dxlink(bed_file),
                  'Contaminants': [dxpy.dxlink(contam_contig)],
                  'Original Reads': [dxpy.dxlink(reads)] }

        print "Running program with", input
        job = self.program.run(input)
        print "Waiting for test_unpaired_with_contam ", job.get_id()
        job.wait_on_done()
        print json.dumps(job.describe()["output"])

        
if __name__ == '__main__':
    unittest.main()
