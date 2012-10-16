

import dxpy
import subprocess
import math
import os

# main function will spawn each analysis as a separate job
# each analysis will then download the mappings and BED gene model and run itself

def run_shell( cmd ):
    print "Running: " + cmd
    subprocess.check_call( cmd, shell=True )

def map_contaminant(Contig, Reads):
    # get ID of our mapper
    try:
        bwa = dxpy.DXApp(dxpy.find_apps(name="bwa").next()['id'])
    except StopIteration:
        raise dxpy.AppError("Unable to find app 'bwa'.  Please install it to enable contaminant mapping")

    map_job = bwa.run({"reads":Reads, "reference": Contig, "discard_unmapped_rows":True, "chunk_size":1000000})

    total_reads = 0
    for r in Reads:
        desc = dxpy.DXGTable(r).describe()
        current_reads = desc['length']
        if 'sequence2' in desc['columns']:
            current_reads *= 2
        total_reads += current_reads

    # launch a job to wait for the mapping and will calculate what % has mapped
    calc_job = dxpy.new_dxjob({"num_reads":total_reads, "mappings":{"job":map_job.get_id(), "field":"mappings"}}, "calc_contam")

    return calc_job.get_id()

@dxpy.entry_point("calc_contam")
def calc_contam(num_reads, mappings):
    percent_mapped = float(dxpy.DXGTable(mappings).describe()['length'])/float(num_reads)

    return {"percent_mapped":percent_mapped}

@dxpy.entry_point("geneBody_coverage")
def geneBody_coverage(BAM_file, BED_file):
    dxpy.download_dxfile(BED_file, "genes.bed")
    dxpy.download_dxfile(BAM_file, "mappings.bam")

    # split mappings into chunks that can be done on a single worker
    # all mappings are loaded into RAM so can only do 10 million at a time
    run_shell(" ".join(["samtools", "view", "mappings.bam", "|", "split", "-l 10000000", "-", "split_map"]))
    run_shell(" ".join(["samtools", "view", "-H", "mappings.bam", ">", "header_only.sam"]))
    files = os.listdir(".")
    jobs = []
    for f in files:
        if f.startswith("split_map"):
            # add header 
            run_shell(" ".join(["cat", "header_only.sam", f, ">", "temp.sam"]))
            # convert to BAM
            run_shell(" ".join(["samtools", "view", "-S", "-b", "temp.sam", ">", "temp.bam"]))
            # upload file
            dxfh = dxpy.upload_local_file("temp.bam")
            # run analysis
            jobs.append(dxpy.new_dxjob({"BAM_file":BAM_file, "BED_file":BED_file}, "run_gbc"))
            
    run_shell( "ls -l" )

    gbc_agg_input = {"sub_reports":[]}
    for j in jobs:
        gbc_agg_input["sub_reports"].append({"job":j.get_id(), "field":"file"})

    agg_job = dxpy.new_dxjob(gbc_agg_input, "gbc_agg").get_id()
    
    return {"results":{"job":agg_job, "field":"cover"}}

@dxpy.entry_point("run_gbc")
def run_gbc(BAM_file, BED_file):
    dxpy.download_dxfile(BED_file, "genes.bed")
    dxpy.download_dxfile(BAM_file, "mappings.bam")

    run_shell( " ".join(["geneBody_coverage.py", "-i mappings.bam", "-r genes.bed", "-o geneBody"]))

    results_id = dxpy.upload_local_file("geneBody.geneBodyCoverage.txt", wait_on_close=True).get_id()
    return {"file":results_id}

@dxpy.entry_point("gbc_agg")
def gbc_agg(sub_reports):
    cover = [0 for n in range(100)]
    total_reads = 0
    for i in range(len(sub_reports)):
        dxpy.download_dxfile(sub_reports[i], str(i)+".txt")
        with open(str(i)+".txt", "r") as fh:
            # remove header info
            fh.readline()
            fh.readline()
            fh.readline()

            for bucket in range(100):
                line = fh.readline()
                count = float(line.split("\t")[1])
                cover[bucket] += count
                total_reads += count

    # normalize by total reads for %
    for i in range(len(cover)):
        cover[i] /= total_reads
    
    return {"cover":cover}

@dxpy.entry_point("inner_distance")
def inner_distance(BAM_file, BED_file):
    dxpy.download_dxfile(BED_file, "genes.bed")
    dxpy.download_dxfile(BAM_file, "mappings.bam")

    run_shell( " ".join(["inner_distance.py", "-i mappings.bam", "-r genes.bed", "-o inner", "-l -303", "-u 5002"]))
    
    results_id = dxpy.upload_local_file("inner.inner_distance_freq.txt", wait_on_close=True).get_id()
    return {"results":results_id}


@dxpy.entry_point("junction_annotation")
def juntion_annotation(BAM_file, BED_file):
    dxpy.download_dxfile(BED_file, "genes.bed")
    dxpy.download_dxfile(BAM_file, "mappings.bam")

    run_shell( " ".join(["junction_annotation.py", "-i mappings.bam", "-r genes.bed", "-o junc"]))

    run_shell( "ls -l" )

    results_id = dxpy.upload_local_file("junc.junction_plot.r", wait_on_close=True).get_id()
    return {"results":results_id}

@dxpy.entry_point("read_distribution")
def read_distribution(BAM_file, BED_file):
    dxpy.download_dxfile(BAM_file, "mappings.bam")
    dxpy.download_dxfile(BED_file, "genes.bed")

    run_shell(" ".join(["read_distribution.py", "-i mappings.bam", "-r genes.bed", ">", "read_dist.txt"]))

    results_id = dxpy.upload_local_file("read_dist.txt", wait_on_close=True).get_id()
    return {"results":results_id}

@dxpy.entry_point("read_duplication")
def read_duplication(BAM_file):
    dxpy.download_dxfile(BAM_file, "mappings.bam")

    run_shell( " ".join(["read_duplication.py", "-i mappings.bam", "-o read_dup"]))
    run_shell( " ".join(["cat", "read_dup.pos.DupRate.xls", "read_dup.seq.DupRate.xls", ">", "read_dup.txt"]))
    results_id = dxpy.upload_local_file("read_dup.txt", wait_on_close=True).get_id()
    return {"results":results_id}

@dxpy.entry_point("generate_report")
def generate_report(geneBody, inner_dist, junc_ann, read_dist, read_dup, mappings, contam, names):
    
    report_details = {}

    # Gene Body Dist
    loc_in_gene = [n for n in range(100)]
    
    report_details['Gene Body Coverage'] = { "Normalized Location in Gene":loc_in_gene,
                                             "% of Reads Covering":geneBody }

    #########################

    # Inner Distance

    if inner_dist != None:

        dxpy.download_dxfile(inner_dist, "inner_dist.txt")

        inner_bucket = []
        inner_num_reads = []
        inner_total_reads = 0
        # if a bucket has less than 0.1% of reads in it then don't include it
        cutoff = 0.001

        with open("inner_dist.inner_distance_freq.txt", "r") as fh:
            line = fh.readline().rstrip("\n")
            while line != "":
                inner_total_reads += int(line.split()[2])
                line = fh.readline().rstrip("\n")

        bucket_cutoff = cutoff * inner_total_reads
        print "Applying cutoff of: "+str(cutoff)+" for inner distance calculation"

        with open("inner_dist.inner_distance_freq.txt", "r") as fh:
            line = fh.readline().rstrip("\n")
            while line != "":
                start, end, num_reads = [int(x) for x in line.split()]
                if num_reads > bucket_cutoff:
                    # store center position of this bucket
                    inner_bucket.append(int(end-((end-start)/2)))
                    inner_num_reads.append(num_reads)

                line = fh.readline().rstrip("\n")

        # find total to normalize
        inner_total_reads = sum(inner_num_reads)
        print "Total reads for inner distance calculation: "+str(inner_total_reads)
        inner_median = None
        running_total = 0
        inner_length_sum = 0
        for i in range(len(inner_bucket)):
            # multiply read length by number of observations for the mean
            inner_length_sum += inner_bucket[i] * inner_num_reads[i]

            # calculate median
            running_total += inner_num_reads[i]
            if running_total >= inner_total_reads / 2 and inner_median == None:
                inner_median = inner_bucket[i]

        inner_mean = inner_length_sum / inner_total_reads
        print "inner distance metrics: "+" ".join([str(inner_length_sum), str(inner_total_reads)])

        # calc standard deviation
        std_sum = 0
        for i in range(len(inner_bucket)):
            std_sum = ((inner_bucket[i] - inner_mean) ** 2) * inner_num_reads[i]

        std_sum /= inner_total_reads
        inner_std = math.sqrt(std_sum)

        report_details['Paired Read Inner Distance'] = {"Inner Distance (bp)": inner_bucket, 
                                                        "Count": inner_num_reads,
                                                        "Mean": inner_mean,
                                                        "Median": inner_median,
                                                        "Standard Deviation": inner_std}

    ############################

    # Junction Annotation

    dxpy.download_dxfile(junc_ann, "junc_ann.r")

    with open("junc_ann.r", "r") as fh:

        line = fh.readline()
        while line != "":
            line = line.rstrip("\n")
            if line.startswith("events"):
                # parse out the % and assign them
                se_pn, se_cn, se_k = [float(n)/100 for n in line[9:-1].split(",")]

            if line.startswith("junction"):
                sj_pn, sj_cn, sj_k = [float(n)/100 for n in line[11:-1].split(",")]

            line = fh.readline()
                
    report_details['Junction Annotation'] = { "Splicing Junctions": {"known": sj_k, "partial novel": sj_pn, "complete novel": sj_cn},
                                              "Splicing Events": {"known": se_k, "partial novel": se_pn, "complete novel": se_cn}}

    ############################

    # read duplication

    dxpy.download_dxfile(read_dup, "read_dup.txt")

    pos_copy = []
    pos_num_reads = []
    pos_total_reads = 0
    seq_copy = []
    seq_num_reads = []
    seq_total_reads = 0

    with open("read_dup.txt", "r") as fh:
        # pull of first header
        line = fh.readline()
        line = fh.readline()
        # read until we hit the stats for sequence based duplication
        while not line.startswith("Occurrence"):
            c, r = [int(n) for n in line.split()]
            pos_copy.append(c)
            pos_num_reads.append(float(r))
            pos_total_reads += r
            line = fh.readline()

        #get next line to start with the data
        line = fh.readline()
        while line != "":
            c, r = [int(n) for n in line.split()]
            seq_copy.append(c)
            seq_num_reads.append(float(r))
            seq_total_reads += r
            line = fh.readline()

    pos_total_reads = float(pos_total_reads)
    seq_total_reads = float(seq_total_reads)

    for i in range(len(pos_num_reads)):
        pos_num_reads[i] /= pos_total_reads

    for i in range(len(seq_num_reads)):
        seq_num_reads[i] /= seq_total_reads

    report_details['Read Duplication'] = {"Position Based":{"Read Occurrences": pos_copy, "% Reads":pos_num_reads},
                                          "Sequence Based":{"Read Occurrences": seq_copy, "% Reads":seq_num_reads}}
    
    ############################

    # read distribution report
    if read_dist != None:
        dxpy.download_dxfile(read_dist, "read_dist.txt")

        report_details['Read Distribution'] = {}

        with open("read_dist.txt", "r") as rd_file:
            report_details['Read Distribution']['Total Reads'] = int(rd_file.readline().split("\t")[1])
            report_details['Read Distribution']['Total Tags'] = int(rd_file.readline().split("\t")[1])
            report_details['Read Distribution']['Total Assigned Tags'] = int(rd_file.readline().split("\t")[1])

            # pull out line of "="s
            rd_file.readline()
            # pull header line
            rd_file.readline()
            line = rd_file.readline()
            while not line.startswith("="):
                fields = line.split("\t")
                report_details['Read Distribution'][fields[0]] = [int(field[1]), int(field[2]), float(field[3])]
                line = rd_file.readline()

    #############################

    # add report of contaminations if calculated

    if contam != None:
        contam_report = []
        for i in range(len(contam)):
            contam_report.append({"Contaminant Name":names[i], "% Reads Mapping":contam[i]})

        report_details['Contamination'] = contam_report
                       
    #############################

    # add link to mappings
    report_details['original_mappings'] = mappings

    report_name = dxpy.DXGTable(mappings).describe()['name'] + " RSeQC report"

    # create report
    report = dxpy.new_dxrecord(name=report_name, details=report_details, types=["Report", "RSeQC"])
    report.close()

    return {"Report": dxpy.dxlink(report.get_id())}

@dxpy.entry_point("main")
def main(**job_inputs):
    output = {}
    reportInput = {}
    
    bed_id = job_inputs["BED file"]
    mappings_id = job_inputs["RNA-Seq Mappings"]["$dnanexus_link"]


    # output mappings as SAM for analysis modules
    run_shell(" ".join(["dx-mappings-to-sam", "--output mappings.sam", mappings_id]))
    run_shell(" ".join(["samtools", "view", "-S", "-b", "mappings.sam", ">", "mappings.bam"]))
    bam_id = dxpy.upload_local_file("mappings.bam", wait_on_close=True).get_id()

    job1 = dxpy.new_dxjob( {'BED_file':bed_id, "BAM_file":dxpy.dxlink(bam_id)}, "geneBody_coverage" )

    # if paired then do inner distance calculation
    if "chr2" in dxpy.DXGTable(mappings_id).get_col_names():
        job2 = dxpy.new_dxjob( {'BED_file':bed_id, "BAM_file":dxpy.dxlink(bam_id)}, "inner_distance" )
    else:
        job2 = None

    job3 = dxpy.new_dxjob( {'BED_file':bed_id, "BAM_file":dxpy.dxlink(bam_id)}, "junction_annotation" )

    job4 = dxpy.new_dxjob( {"BAM_file":dxpy.dxlink(bam_id)}, "read_duplication" )

    # implement this one when we can request a large RAM instance - requires 19GB for human genome
    #job5 = dxpy.new_dxjob( {'BED_file':bed_id, "BAM_file":dxpy.dxlink(bam_id)}, "read_distribution", 
    #                       {"systemRequirements": {"instanceType":"dx_m2.2xlarge"}} )

    # get contaminant mapping started if we're doing it:
    if "Contaminants" in job_inputs:
        name_input = []
        contam_input = []

        #spawn mappings job for each ContigSet
        for contaminant in job_inputs['Contaminants']:
            calc_job = map_contaminant(Reads=job_inputs['Original Reads'], Contig=contaminant)

            name_input.append(dxpy.DXRecord(contaminant).describe()['name'])
            contam_input.append({"job":calc_job, "field":"percent_mapped"})
    
        reportInput['contam'] = contam_input
        reportInput['names'] = name_input
    else:
        reportInput['contam'] = None
        reportInput['names'] = None


    reportInput['geneBody'] = {"job":job1.get_id(), "field":"results"}
    if job2 != None:
        reportInput['inner_dist'] = {"job":job2.get_id(), "field":"results"}
    else:
        reportInput['inner_dist'] = None
    reportInput['junc_ann'] = {"job":job3.get_id(), "field":"results"}
    reportInput['read_dup'] = {"job":job4.get_id(), "field":"results"}
    #reportInput['read_dist'] = {"job":job5.get_id(), "field":"results"}
    reportInput['read_dist'] = None
    reportInput['mappings'] = job_inputs["RNA-Seq Mappings"]


    reportJob = dxpy.new_dxjob( reportInput, "generate_report" )
    

    output['Report'] = {"job":reportJob.get_id(), "field": "Report"}
    
    return output
