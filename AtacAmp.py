import argparse
import os
import re
import subprocess
import sys
import time
from multiprocessing import Pool

import interval
import psutil
import pysam

import calculate_cnv


def parse_args():
    parser = argparse.ArgumentParser(description='use atac data find ecdna')
    parser.add_argument('--bam', type=str, dest='bam', action='store', help='input the bam file')
    parser.add_argument('--name', '-n', type=str, dest='name', action='store', help='prefix of output files')
    parser.add_argument('--isize_value', '-i', type=int, dest='isize', action='store', help='judge a pair of reads whether is discordant')
    parser.add_argument('--interval_size', '-s', type=int, dest='interval', action='store', help='size of interval when compute breakpoint nearby coverage')
    parser.add_argument('--mapq', '-q', type=int, dest='maqp', action='store', help='reads maqp threshold')
    parser.add_argument('--mode', '-m', type=int, dest='mode', action='store', choices=[0, 1, 2], help='choose the analysis mode,0/1/2')
    parser.add_argument('--discbk', '-d', type=str, dest='discbk', action='store', help='if you choose mode 1,you need to input discbk file')
    parser.add_argument('--type', type=str, dest='lib', action='store', choices=['sc', 'bulk'], help='choose library:sc/bulk')
    parser.add_argument('--gtf', type=str, dest='gtf', action='store', help='gtf file')
    parser.add_argument('--threads', type=int, dest='threads', action='store', default=1, help='threads')

    args = parser.parse_args()

    if args.mode == 2:
        parser.error('mode 2 is not implemented in current version')
    if not args.bam:
        parser.error('--bam is required')
    if args.mode == 0 and args.isize is None:
        parser.error('--isize_value/-i is required in mode 0')
    if args.mode == 1 and not args.discbk:
        parser.error('--discbk/-d is required in mode 1')
    if args.interval is None:
        parser.error('--interval_size/-s is required')
    if not args.gtf:
        parser.error('--gtf is required')
    if not args.lib:
        parser.error('--type is required and must be sc or bulk')

    if args.name is None:
        args.name = os.path.basename(args.bam).rstrip('.bam')

    return args


def calculate_average_depth_for_contig(args):
    bam_file, contig = args
    bam = pysam.AlignmentFile(bam_file, 'rb')
    total_depth = 0
    total_positions = 0
    for pileup_column in bam.pileup(contig):
        total_depth += pileup_column.n
        total_positions += 1
    bam.close()
    return total_depth, total_positions


def calculate_average_depth(bam_file, num_threads=1):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    contigs = bam.references
    bam.close()

    with Pool(num_threads) as pool:
        results = pool.map(calculate_average_depth_for_contig, [(bam_file, contig) for contig in contigs])

    total_depth = 0
    total_positions = 0
    for depth, positions in results:
        total_depth += depth
        total_positions += positions

    return total_depth / total_positions if total_positions else 0


def open_alignment_file(bam_file):
    if bam_file.endswith('.bam'):
        return pysam.AlignmentFile(bam_file, 'rb')
    if bam_file.endswith('.sam'):
        return pysam.AlignmentFile(bam_file, 'r')
    raise ValueError('please input sam or bam file')


def find_special_reads(bam_file, pref, isize, mapq_filter=0):
    f_out_bam_split = pref + '.split.bam'
    f_out_bam_discordant = pref + '.discordant.bam'

    input_bam = open_alignment_file(bam_file)
    out_bam_split = pysam.AlignmentFile(f_out_bam_split, 'wb', template=input_bam)
    out_bam_discordant = pysam.AlignmentFile(f_out_bam_discordant, 'wb', template=input_bam)

    line_n = 0
    s_t = time.time()

    for reads_line in input_bam:
        line_n += 1
        if line_n % 2000000 == 0:
            print(line_n)
            print(time.time() - s_t)
            print('memory：%.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024))

        if (reads_line.reference_name in ['chrM', 'MT']) or (reads_line.next_reference_name in ['chrM', 'MT']):
            continue

        if mapq_filter and reads_line.mapq < mapq_filter:
            continue

        if reads_line.cigarstring is None:
            continue

        if ('S' in reads_line.cigarstring) or ('H' in reads_line.cigarstring):
            out_bam_split.write(reads_line)

        if (abs(reads_line.isize) > isize) or (reads_line.next_reference_name != reads_line.reference_name):
            out_bam_discordant.write(reads_line)

    input_bam.close()
    out_bam_split.close()
    out_bam_discordant.close()

    return f_out_bam_split, f_out_bam_discordant


def _cigar_breakpoint_offset(cigar_prefix):
    matches = re.findall(r'(\d+)([IDM])', cigar_prefix)
    sum_insertions = 0
    sum_deletions = 0
    sum_matches = 0
    for length, operation in matches:
        length = int(length)
        if operation == 'I':
            sum_insertions += length
        elif operation == 'D':
            sum_deletions += length
        elif operation == 'M':
            sum_matches += length
    return sum_matches + sum_deletions - sum_insertions


def process_split_reads(split_bam):
    bam = pysam.AlignmentFile(split_bam, 'rb')

    outfile = split_bam.split('/')[-1].rstrip('.bam') + '.split_bk'
    f_out = open(outfile, 'w')

    breakpoint_dir = {}

    s_t = time.time()
    count1 = 0
    for split_line in bam:
        count1 += 1
        if count1 % 10000 == 0:
            print(count1)
            print(time.time() - s_t)

        if not split_line.has_tag('SA'):
            continue

        read_cigar_prefix = split_line.cigarstring.split('M')[0]
        if not (('H' in read_cigar_prefix) or ('S' in read_cigar_prefix)):
            breakpoint_start = split_line.pos + int(_cigar_breakpoint_offset(read_cigar_prefix))
            orient1 = '+'
        else:
            breakpoint_start = split_line.pos
            orient1 = '-'

        sa_info = split_line.get_tag('SA').split(',')
        sa_chr = sa_info[0]
        sa_start = sa_info[1]
        sa_cigar = sa_info[3]

        sa_cigar_prefix = sa_cigar.split('M')[0]
        if not (('H' in sa_cigar_prefix) or ('S' in sa_cigar_prefix)):
            breakpoint_end = int(sa_start) + int(_cigar_breakpoint_offset(sa_cigar_prefix))
            orient2 = '+'
        else:
            breakpoint_end = int(sa_start)
            orient2 = '-'

        read1 = split_line.reference_name + '\t' + str(breakpoint_start) + '\t' + orient1
        read2 = sa_chr + '\t' + str(breakpoint_end) + '\t' + orient2

        if (split_line.reference_name > sa_chr) or (breakpoint_end <= breakpoint_start):
            read1, read2 = read2, read1
        breakpoint_info = read1 + '\t' + read2

        if breakpoint_info in breakpoint_dir:
            breakpoint_dir[breakpoint_info].append(split_line)
        else:
            breakpoint_dir[breakpoint_info] = [split_line]

    for key, value in breakpoint_dir.items():
        f_out.write(key + '\t' + str(len(value)) + '\n')

    f_out.close()
    bam.close()


class BreakpointObj:
    def __init__(self, chr_l, l_pos, chr_r, r_pos, qname, flag, support=1):
        self.chr_l = chr_l
        self.chr_r = chr_r
        self.l_pos = l_pos
        self.r_pos = r_pos
        self.support_reads = [qname + '\t' + str(flag)]
        self.support = support

    def modify(self, new_l, new_r, new_qname, new_flag):
        self.support += 1
        self.support_reads.append(new_qname + '\t' + str(new_flag))
        self.l_pos = new_l
        self.r_pos = new_r

    def integrate(self):
        return self.chr_l + '\t' + str(self.l_pos) + '\t' + self.chr_r + '\t' + str(self.r_pos) + '\t' + str(self.support) + '\t' + '\t'.join(self.support_reads)


def process_discordant_reads(disc_bam_file, f_long=500):
    disc_bam = open_alignment_file(disc_bam_file)
    disc_outfile = os.path.basename(disc_bam_file.rsplit('.', 1)[0])

    f_disc_breakpoint = open(disc_outfile + '.disc_bk', 'w')
    breakpoint_dir = {ref: [] for ref in disc_bam.references}

    count2 = 0
    s_t = time.time()

    for disc_line in disc_bam:
        count2 += 1
        if count2 % 10000 == 0:
            print(time.time() - s_t)
            print(count2)

        if disc_line.reference_name == disc_line.next_reference_name:
            bk_start = min(disc_line.pos, disc_line.next_reference_start)
            bk_end = max(disc_line.pos, disc_line.next_reference_start)

            temp_list = [
                interval.Interval((bk_start + 1) - f_long, (bk_start + 1) + f_long),
                interval.Interval((bk_end + 1) - f_long, (bk_end + 1) + f_long),
                disc_line.reference_name,
                disc_line.qname,
                disc_line.flag,
            ]
            disc_bk_line = BreakpointObj(temp_list[2], temp_list[0], temp_list[2], temp_list[1], temp_list[3], temp_list[4])

            i = 1
            for list_line in breakpoint_dir[disc_line.reference_name]:
                if list_line.l_pos.overlaps(temp_list[0]) and list_line.r_pos.overlaps(temp_list[1]):
                    list_line.modify(list_line.l_pos.join(temp_list[0]), list_line.r_pos.join(temp_list[1]), temp_list[3], temp_list[4])
                    i = 0
                    break

            if i == 1:
                breakpoint_dir[disc_line.reference_name].append(disc_bk_line)

        else:
            breakpoint_r = interval.Interval((disc_line.next_reference_start + 1 - f_long), (disc_line.next_reference_start + 1 + f_long))
            breakpoint_l = interval.Interval((disc_line.pos + 1 - f_long), (disc_line.pos + 1 + f_long))

            if disc_line.reference_name <= disc_line.next_reference_name:
                temp_chr_l = disc_line.reference_name
                temp_chr_r = disc_line.next_reference_name
            else:
                temp_chr_l = disc_line.next_reference_name
                temp_chr_r = disc_line.reference_name
                breakpoint_l, breakpoint_r = breakpoint_r, breakpoint_l

            temp_key = temp_chr_l + '\t' + temp_chr_r

            if temp_key in breakpoint_dir:
                index_temp1 = 1
                for bk_line in breakpoint_dir[temp_key]:
                    if bk_line.l_pos.overlaps(breakpoint_l) and bk_line.r_pos.overlaps(breakpoint_r):
                        bk_line.modify(bk_line.l_pos.join(breakpoint_l), bk_line.r_pos.join(breakpoint_r), disc_line.qname, disc_line.flag)
                        index_temp1 = 0
                        break
                if index_temp1 == 1:
                    breakpoint_dir[temp_key].append(BreakpointObj(temp_chr_l, breakpoint_l, temp_chr_r, breakpoint_r, disc_line.qname, disc_line.flag))
            else:
                breakpoint_dir[temp_key] = [BreakpointObj(temp_chr_l, breakpoint_l, temp_chr_r, breakpoint_r, disc_line.qname, disc_line.flag)]

    for value in breakpoint_dir.values():
        for line in value:
            f_disc_breakpoint.write(line.integrate() + '\n')

    disc_bam.close()
    f_disc_breakpoint.close()
    return disc_outfile + '.disc_bk'


def run_command(cmd):
    print(' '.join(cmd))
    subprocess.run(cmd, check=True)


def run_post_processing(args, bkamplicon_f, disc_bam_f):
    argv2 = args.name + '.bkline_interval'
    argv3 = args.name + '.bkline_dif_interval'

    run_command([sys.executable, os.path.join(sys.path[0], 'bkgraph.py'), bkamplicon_f, argv2, argv3])
    print('bkgraph finished')

    run_command(['samtools', 'index', disc_bam_f])
    run_command([
        sys.executable,
        os.path.join(sys.path[0], '20230301_graph.py'),
        argv3,
        args.name + '.result',
        args.lib,
        args.gtf,
        str(args.threads),
    ])


def run_mode0(args):
    split_bam_f, disc_bam_f = find_special_reads(args.bam, args.name, args.isize, args.maqp or 0)
    print('find special reads finished')

    process_split_reads(split_bam_f)
    print('process split_reads finished')

    discbk_f = process_discordant_reads(disc_bam_f)
    print('process_disc_reads_finished')

    acov = calculate_average_depth(args.bam, args.threads)
    bkamplicon_f = calculate_cnv.find_amp(args.bam, discbk_f, interval=args.interval, cov=acov).process_coverage()
    print('amplicon find finished')

    run_post_processing(args, bkamplicon_f, disc_bam_f)


def run_mode1(args):
    s_t = time.time()
    acov = calculate_average_depth(args.bam, args.threads)

    bkamplicon_f = calculate_cnv.find_amp(args.bam, args.discbk, interval=args.interval, cov=acov).process_coverage()
    print('amplicon find finished')
    print(time.time() - s_t)

    disc_bam_f = args.discbk.replace('disc_bk', 'bam')
    run_post_processing(args, bkamplicon_f, disc_bam_f)


def main():
    args = parse_args()

    if args.mode == 0:
        run_mode0(args)
    elif args.mode == 1:
        run_mode1(args)


if __name__ == '__main__':
    main()
