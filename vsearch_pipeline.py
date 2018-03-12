#!/usr/bin/python
import sys, os, argparse, string
import glob
from Bio import SeqIO
from subprocess import call, Popen, PIPE

# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description='')
requiredArguments = parser.add_argument_group('required arguments')
requiredArguments.add_argument('-i', '--input', metavar='input zipfile', dest='inzip', type=str,
                               help='Inputfile in zip format', required=True)
requiredArguments.add_argument('-o', '--output_folder', metavar='output', dest='out', type=str,
                               help='output in zip format', required=True, nargs='?', default="")
requiredArguments.add_argument('-t', '--input_type', metavar='FASTQ or FASTA input', dest='input_type', type=str,
                               help='Sets the input type, FASTQ or FASTA', required=True)
requiredArguments.add_argument('-cluster_id', '--cluster_id', metavar='Minimal cluster identity percentage', dest='clusterid', type=str,
                               help='Minimal cluster identity percentage', required=True)
requiredArguments.add_argument('-cluster_size', '--cluster_size', metavar='Minimal cluster size', dest='clustersize', type=str,
                               help='Minimal cluster size', required=True)
args = parser.parse_args()

def admin_log(out=None, error=None, function=""):
    with open(args.out+"/adminlog.log", 'a') as adminlogfile:
        seperation = 60 * "="
        if out:
            adminlogfile.write("out "+ function + " \n" + seperation + "\n" + out + "\n\n")
        if error:
            adminlogfile.write("error " + function + "\n" + seperation + "\n" + error + "\n\n")

def check_if_fasta(file):
    with open(file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def format_check():
    files = [os.path.basename(x) for x in sorted(glob.glob(args.out+"/inputfiles/*"))]
    for x in files:
        if args.input_type == "FASTQ":
            if os.path.splitext(x)[1].lower() == ".fastq" or os.path.splitext(x)[1] == ".fq":
                fastafile = os.path.splitext(x)[0].translate((string.maketrans("-. " , "___")))+".fa"
                error = Popen(["awk '{if(NR%4==1) {printf(\">%s\\n\",substr($0,2));} else if(NR%4==2) print;}' "+ args.out+"/inputfiles/" + x + " > "+args.out+"/fastafiles/" + fastafile], stdout=PIPE, stderr=PIPE, shell=True).communicate()[1].strip()
                admin_log(error=error, function="extension_check")
                call(["sed 's/>/>" + fastafile[:-3] + "./' " + args.out + "/fastafiles/"+fastafile+" >>" + args.out+ "/workfiles/combined_raw.fa"], shell=True)
                out, error = Popen(["vsearch", "--derep_fulllength", args.out+ "/fastafiles/"+fastafile, "--sizeout", "--fasta_width", "0", "--relabel", fastafile[:-3]+".","--output", args.out+ "/fastafiles/"+fastafile[:-3]+".derep.fa","--sizeout"], stdout=PIPE, stderr=PIPE).communicate()
                admin_log(out=out, error=error, function="dereplicate per sample")
                call(["rm", args.out+ "/fastafiles/"+fastafile])
            else:
                admin_log(error=x+"\nWrong extension, no fastq file (.fastq, .fq) file will be ignored", function="format_check")
        else:
            if check_if_fasta(args.out+ "/inputfiles/" + x):
                fastafile = os.path.splitext(x)[0].translate((string.maketrans("-. ", "___"))) + ".fa"
                call(["mv", args.out+ "/inputfiles/" + x, args.out+ "/fastafiles/" + fastafile])
                call(["sed 's/>/>" + fastafile[:-3] + "./' "+args.out+ "/fastafiles/" + fastafile + " >> " +args.out+ "/workfiles/combined_raw.fa"], shell=True)
                out, error = Popen(["vsearch", "--derep_fulllength", args.out+ "/fastafiles/"+fastafile, "--sizeout", "--fasta_width", "0", "--relabel", fastafile[:-3]+".","--output", args.out+ "/fastafiles/"+fastafile[:-3]+".derep.fa","--sizeout"], stdout=PIPE, stderr=PIPE).communicate()
                admin_log(out=out, error=error, function="dereplicate per sample")
                call(["rm", args.out+ "/fastafiles/" + fastafile])
            else:
                admin_log(error="This is not a fasta file, file will be ignored: " + x, function="extension_check")
    Popen(["rm", "-rf", args.out+ "/inputfiles"], stdout=PIPE, stderr=PIPE)

def vsearch_derep_fulllength():
    call(["cat " +args.out+ "/fastafiles/*.derep.fa > " +args.out+ "/workfiles/combined.fa"], shell=True)
    out, error = Popen(["vsearch", "--derep_fulllength", args.out+ "/workfiles/combined.fa", "--output", args.out+ "/workfiles/uniques.fa", "--sizein","--sizeout","--fasta_width", "0"], stdout=PIPE, stderr=PIPE).communicate()
    admin_log(out=out, error=error, function="dereplicate merged samples")

def chimera_check():
    out, error = Popen(["vsearch", "--uchime_denovo", args.out+ "/workfiles/uniques.fa", "--sizein", "--fasta_width", "0", "--nonchimeras", args.out+ "/workfiles/non_chimera.fa"], stdout=PIPE, stderr=PIPE).communicate()
    admin_log(out=out, error=error, function="chimera_check")

def cluster():
    out, error = Popen(["vsearch", "--cluster_size", args.out+ "/workfiles/non_chimera.fa", "--id", args.clusterid, "--sizein", "--sizeout", "--fasta_width", "0", "--relabel", "OTU_", "--centroids", args.out+ "/workfiles/otu_sequences.fa"], stdout=PIPE,stderr=PIPE).communicate()
    admin_log(out=out, error=error, function="clustering")

def make_otu_table():
    out, error = Popen(["vsearch", "--usearch_global", args.out+ "/workfiles/combined_raw.fa", "--db", args.out+ "/workfiles/otu_sequences.fa", "--id", "0.98", "--otutabout", args.out+ "/workfiles/otutab.txt"], stdout=PIPE, stderr=PIPE).communicate()
    admin_log(out=out, error=error, function="make_otu_table")

def remove_clusters():
    otuCount = 0
    singletonCount = 0
    with open(args.out+ "/workfiles/filtered_otu.fa", "a") as otuFile:
        for record in SeqIO.parse(args.out + "/workfiles/otu_sequences.fa", "fasta"):
            size = str(record.description).split("size=")[1][:-1]
            if int(size) > int(args.clustersize):
                otuFile.write(">"+str(record.description)+"\n"+str(record.seq)+"\n")
                otuCount += 1
            else:
                singletonCount += 1
    #write to log
    out = "Minimal otu size:" + str(args.clustersize) + " \nOtu's:" + str(otuCount) + "\nSingletons:" + str(singletonCount)
    admin_log(out=out, function="cluster size filtering")

def main():
    args.out = str(args.out).strip()
    #make folders to store files
    call(["mkdir", args.out+"/inputfiles", args.out+"/fastafiles", args.out+"/workfiles"])
    #upzip the input file, write the files to the inputfiles folder
    zip_out, zip_error = Popen(["unzip", args.inzip, "-d", args.out+"/inputfiles"], stdout=PIPE,stderr=PIPE).communicate()
    #write the streams of the zip command to the log file
    admin_log(out=zip_out, error=zip_error, function="unzip")
    #check the input format, this can be fasta or fastq and dereplicate each sample
    format_check()
    #dereplication of all samples merged
    vsearch_derep_fulllength()
    #chimera check
    chimera_check()
    #cluster
    cluster()
    #remove cluster based on size treshold
    remove_clusters()
    #make otu table
    make_otu_table()


if __name__ == '__main__':
    main()