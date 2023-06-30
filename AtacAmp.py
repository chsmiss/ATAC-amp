import argparse
import pysam
import time
#import multiprocessing
import psutil
import os
import sys
#import concurrent.futures
import interval
import calculate_cnv
import subprocess
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='use atac data find ecdna')
parser.add_argument('--bam',type=str,dest='bam',action='store',help = 'input the bam file')
parser.add_argument('--name','-n',type=str,dest='name',action='store',help = 'prefix of output files ')
parser.add_argument('--isize_value','-i',type=int,dest='isize',action='store',help = 'judge a pair of reads whether is discordant')
parser.add_argument('--interval_size','-s',type=int,dest='interval',
                    action='store',help = 'size of interval when compute breakpoint nearby coverage')
parser.add_argument('--mapq','-q',type=int,dest='maqp',action='store',help = 'reads maqp threshold')
parser.add_argument('--mode','-m',type=int,dest='mode',action='store',choices=[0, 1, 2],help = 'choose the analysis mode,0/1/2')
parser.add_argument('--discbk','-d',type=str,dest='discbk',action='store',help = 'if you choose mode 1,you need to input discbk file')
parser.add_argument('--type',type=str,dest='lib',action='store',help = 'choose library:sc/bulk')
parser.add_argument('--gtf',type=str,dest='gtf',action='store',help = 'gtf file')
parser.add_argument('--threads',type=int,dest='threads',action='store',help = 'threads')
args = parser.parse_args()

if args.name == None:
    args.name = args.bam.split('/')[-1].rstrip('.bam')

f_input_bam = args.bam
isize_value = args.isize
f_prefix = args.name
extent_interval_size = args.interval


def calculate_average_depth_for_contig(args):
    bam_file, contig = args
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_depth = 0
    total_positions = 0
    for pileup_column in bam.pileup(contig):
        total_depth += pileup_column.n
        total_positions += 1
    bam.close()
    return total_depth, total_positions

def calculate_average_depth(bam_file, num_threads=1):
    bam = pysam.AlignmentFile(bam_file, "rb")
    contigs = bam.references
    bam.close()
    pool = Pool(num_threads)
    results = pool.map(calculate_average_depth_for_contig, [(bam_file, contig) for contig in contigs])
    total_depth = 0
    total_positions = 0
    for depth, positions in results:
        total_depth += depth
        total_positions += positions
    average_depth = total_depth / total_positions
    return average_depth

    bam.close()
    average_depth = total_depth / total_positions
    return average_depth

#def find_special_reads(bam_file,pref,isize,mapq_filter=0):
def find_special_reads(bam_file, pref, isize,mapq_filter=0):

    f_out_bam_split = pref + '.split.bam'
    f_out_bam_discordant = pref+ '.discordant.bam'

    if f_input_bam.endswith('.bam'):
        input_bam = pysam.AlignmentFile(bam_file, 'rb')
    elif f_input_bam.endswith('.sam'):
        input_bam = pysam.AlignmentFile(bam_file, 'r')
    else:
        print('please input sam or bam file')

    out_bam_split = pysam.AlignmentFile(f_out_bam_split, 'wb', template=input_bam)
    out_bam_discordant = pysam.AlignmentFile(f_out_bam_discordant, 'wb', template=input_bam)

    line_n = 0
    s_t = time.time()

    for reads_line in input_bam:
        line_n = line_n + 1
        if line_n%2000000 == 0:
            print(line_n)
            e_t = time.time()
            print(e_t - s_t)
            print('membrane：%.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024))

        if (reads_line.reference_name in ['chrM', 'MT']) or (reads_line.next_reference_name in ['chrM', 'MT']):
            continue

        if mapq_filter:
            if reads_line.mapq < mapq_filter:
                continue

        if reads_line.cigarstring == None:
            continue

        if ('S' in reads_line.cigarstring) or ('H' in reads_line.cigarstring):
            out_bam_split.write(reads_line)


        if (abs(reads_line.isize) > isize_value) or (reads_line.next_reference_name != reads_line.reference_name):
            out_bam_discordant.write(reads_line)

    input_bam.close()
    out_bam_split.close()
    out_bam_discordant.close()

    return (f_out_bam_split,f_out_bam_discordant)


def process_split_reads(split_bam):
    bam = pysam.AlignmentFile(split_bam, 'rb')

    outfile = split_bam.split('/')[-1].rstrip('.bam')+'.split_bk'
    f_out = open(outfile, 'w')

    breakpoint_dir = {}

    s_t = time.time()
    count1 = 0
    for split_line in bam:

        count1=count1+1
        if count1 % 10000 == 0:
            print(count1)
            print(time.time()-s_t)

        if split_line.has_tag('SA'):

            if not (('H' in split_line.cigarstring.split('M')[0]) or ('S' in split_line.cigarstring.split('M')[0])):
                breakpoint_start = split_line.pos + int(split_line.cigarstring.split('M')[0])
                orient1 = '+'
            else:
                breakpoint_start = split_line.pos
                orient1 = '-'

            sa_info = split_line.get_tag('SA').split(',')
            sa_chr = sa_info[0]
            sa_start = sa_info[1]
#            sa_chain = sa_info[2]
            sa_cigar = sa_info[3]


            sa_cigar_info = sa_cigar.split('M')
            if not (('H' in sa_cigar_info[0]) or ('S' in sa_cigar_info[0])):
                breakpoint_end = int(sa_start) + int(sa_cigar_info[0])
                orient2 = '+'
            else:
                breakpoint_end = int(sa_start)
                orient2 = '-'
            read1 = split_line.reference_name + '\t' + str(breakpoint_start)+'\t'+orient1
            read2 = sa_chr + '\t' + str(breakpoint_end) +'\t'+ orient2

            if (split_line.reference_name > sa_chr) or (breakpoint_end <= breakpoint_start):
       #         breakpoint_info = sa_chr + '\t' + str(breakpoint_end) + '\t' + split_line.reference_name + '\t' + str(breakpoint_start)
                read1,read2 = read2,read1
            breakpoint_info = read1 + '\t' + read2


            if breakpoint_info in breakpoint_dir:
                breakpoint_dir[breakpoint_info].append(split_line)
            else:
                breakpoint_dir[breakpoint_info] = [split_line]


    for key, value in breakpoint_dir.items():
      #  value_line = ''
        f_out.write(key + '\t' + str(len(value)) + '\n')

  #  return (breakpoint_dir)


class BreakpointObj:
    def __init__(self, chr_l, l_pos, chr_r, r_pos, qname, flag, support=1):
        self.chr_l = chr_l
        self.chr_r = chr_r
        self.l_pos = l_pos
        self.r_pos = r_pos
        self.support_reads = []
        self.support_reads.append(qname + '\t' + str(flag))
        self.support = support

    def modify(self, new_l, new_r, new_qname, new_flag):
        self.support = self.support + 1
        self.support_reads.append(new_qname + '\t' + str(new_flag))
        self.l_pos = new_l
        self.r_pos = new_r

    def integrate(self):
        #		break_line = self.chr_l+str(self.l_pos)+self.chr_r+str(self.r_pos)+'\t'+'\t'.join(';'.join(self.suport_reads))
        return self.chr_l + '\t' + str(self.l_pos) + '\t' + self.chr_r + '\t' + str(self.r_pos) + '\t' + str(
            self.support) + '\t' + '\t'.join(self.support_reads)
#		return(break_line)

def sort_interval_list(interval_list):
# for element in interval_list:
#     element.l_pos
        pass


def process_discordant_reads(disc_bam_file,f_long=500):
    if disc_bam_file.endswith('.bam'):
        disc_bam = pysam.AlignmentFile(disc_bam_file, 'rb')
        disc_outfile = disc_bam_file.strip('.bam').split('/')[-1]
    if disc_bam_file.endswith('.sam'):
        disc_bam = pysam.AlignmentFile(disc_bam_file, 'r')
        disc_outfile = disc_bam_file.strip('.sam').split('/')[-1]

    f_disc_breakpoint = open(disc_outfile+'.disc_bk', 'w')

    breakpoint_dir = {}
  #  bk_dir = []


    count2=0
    s_t = time.time()
    for ref in disc_bam.references:
        breakpoint_dir[ref] = []

    for disc_line in disc_bam:
        count2=count2+1
        if count2 % 10000 == 0:
            print(time.time()-s_t)
            print(count2)

        #read1和read2比对到同一染色体上
        if disc_line.reference_name == disc_line.next_reference_name:
            bk_start = min(disc_line.pos, disc_line.next_reference_start)
            bk_end = max(disc_line.pos, disc_line.next_reference_start)


            temp_list = [interval.Interval((bk_start+1)-f_long,(bk_start+1)+f_long),interval.Interval((bk_end+1)-f_long,(bk_end+1)+f_long),disc_line.reference_name,disc_line.qname,disc_line.flag]
            disc_bk_line = BreakpointObj(temp_list[2],temp_list[0],temp_list[2],temp_list[1],temp_list[3],temp_list[4])


            if disc_line.reference_name in breakpoint_dir.keys():
                i = 1
        
   #             for index, list_line in enumerate(breakpoint_dir[disc_line.reference_name]):

                for list_line in breakpoint_dir[disc_line.reference_name]:
                    if list_line.l_pos.overlaps(temp_list[0]) and list_line.r_pos.overlaps(temp_list[1]):
                        list_line.modify(list_line.l_pos.join(temp_list[0]),list_line.r_pos.join(temp_list[1]),temp_list[3],temp_list[4])
                        i = 0
                        break
     
                if i == 1:
                    breakpoint_dir[disc_line.reference_name].append(disc_bk_line)


            # else:
            #     breakpoint_dir[disc_line.reference_name] = []
            #     breakpoint_dir[disc_line.reference_name].append(disc_bk_line)


        else:
            breakpoint_r = interval.Interval((disc_line.next_reference_start+1-f_long), (disc_line.next_reference_start+1+f_long))
            breakpoint_l = interval.Interval((disc_line.pos+1-f_long), (disc_line.pos+1+f_long))

            if disc_line.reference_name > disc_line.next_reference_name:
                temp_chr_l,temp_chr_r = disc_line.next_reference_name,disc_line.reference_name
                breakpoint_l,breakpoint_r = breakpoint_r,breakpoint_l

            temp_chr_r, temp_chr_l = disc_line.next_reference_name, disc_line.reference_name
            temp_key = temp_chr_l+'\t'+temp_chr_r

            if temp_key in breakpoint_dir:
                index_temp1 = 1
                for bk_line in breakpoint_dir[temp_key]:
                    if bk_line.l_pos.overlaps(breakpoint_l) and bk_line.r_pos.overlaps(breakpoint_r):
                        bk_line.modify(bk_line.l_pos.join(breakpoint_l), bk_line.r_pos.join(breakpoint_r), disc_line.qname, disc_line.flag)
                        index_temp1 = 0
                        break
                if index_temp1 == 1:
                    breakpoint_dir[temp_key].append(BreakpointObj(disc_line.next_reference_name, breakpoint_l, disc_line.reference_name,breakpoint_r, disc_line.qname, disc_line.flag))


            else:
                breakpoint_dir[temp_key]=[]
                breakpoint_dir[temp_key].append(BreakpointObj(disc_line.next_reference_name, breakpoint_l, disc_line.reference_name, breakpoint_r, disc_line.qname, disc_line.flag))


    for key, value in breakpoint_dir.items():
        for line in value:
            f_disc_breakpoint.write(line.integrate() + '\n')

    f_disc_breakpoint.close()
  #  return (breakpoint_dir,bk_dir)
    disc_file_name = disc_outfile+'.disc_bk'
    return (disc_file_name)


#test_cnv.process_coverage(f_input_bam)
if args.mode == 0:
    split_bam_f,disc_bam_f = find_special_reads(f_input_bam,f_prefix,isize_value)
    print('find special reads finished')
    split_bk_dir = process_split_reads(split_bam_f)
    process_split_reads(split_bam_f)
    print('process split_reads finished')
    discbk_f = process_discordant_reads(disc_bam_f)
    print('process_disc_reads_finished')
    acov = calculate_average_depth(f_input_bam, args.threads)

    bkamplicon_f = calculate_cnv.find_amp(f_input_bam, discbk_f,
                                          interval=extent_interval_size, cov=
                                          acov).process_coverage()
    print('amplicon find finished')
    argv2 = args.name+'.bkline_interval'
    argv3 = args.name+'.bkline_dif_interval'
    argv_line = 'python '+sys.path[0]+'/bkgraph.py '+bkamplicon_f+' '+argv2+' '+argv3
    subprocess.run(argv_line,shell=True)
    print('bkgraph finished')

    index_line = 'samtools index '+disc_bam_f
    subprocess.run(index_line,shell=True)
    if args.lib == 'sc':
        argv_line2 = 'python '+sys.path[0]+'/20230301_graph.py '+argv3+' '+args.name+'.result '+'sc '+args.gtf + ' '+str(args.threads)
    if args.lib == 'bulk':
        argv_line2 = 'python '+sys.path[0]+'/20230301_graph.py ' + argv3 + ' ' + args.name + '.result ' + 'bulk ' + args.gtf + ' '+str(args.threads)
    print(argv_line2)
    subprocess.run(argv_line2,shell=True)

if args.mode == 1:
    discbk_f = args.discbk
    s_t = time.time()
    acov = calculate_average_depth(f_input_bam, args.threads)

    bkamplicon_f = calculate_cnv.find_amp(f_input_bam, discbk_f,
                              interval=extent_interval_size,cov =
                                          acov).process_coverage()
    print('amplicon find finished')
    print(time.time() - s_t)
    argv2 = args.name + '.bkline_interval'
    argv3 = args.name + '.bkline_dif_interval'
    argv_line = 'python '+sys.path[0]+'/bkgraph.py ' + bkamplicon_f + ' ' + argv2 + ' ' + argv3
    subprocess.run(argv_line, shell=True)

    #20230501 add
    print('bkgraph finished')
    disc_bam_f = discbk_f.strip('disc_bk')+'bam'
    index_line = 'samtools index '+disc_bam_f
    subprocess.run(index_line,shell=True)
    if args.lib == 'sc':
        argv_line2 = 'python '+sys.path[0]+'/20230301_graph.py '+argv3+' ' \
                                                                       ''+args.name+'.result '+'sc '+args.gtf+' '+str(args.threads)
    if args.lib == 'bulk':
        argv_line2 = 'python '+sys.path[0]+'/20230301_graph.py ' + argv3 + ' ' \
                                                                           ''\
                     + args.name + '.result ' + 'bulk ' + args.gtf+' ' \
                                                                   ''+str(args.threads)
    print(argv_line2)
    subprocess.run(argv_line2,shell=True)


