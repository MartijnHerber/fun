#!/usr/bin/env python3
# (c) 2022 Martijn Herber / Stichting Hanzehogeschool Groningen
# Licensed under GPL 3.0; see file LICENSE for details.

def simple_pure_main(bedfile_name, pileup_name):
    from collections import defaultdict, OrderedDict
    bed_fp = open(bedfile_name, "r")
    pileup_fp = open(pileup_name, "r")

    chrms = defaultdict(OrderedDict)
    genes = {}
    lowcov = []
    
    for line in bed_fp:
        chrm, start, stop, gene = line.split(sep="\t")
        chrms["chr"+chrm].update({str(i) : gene for i in range(int(start), int(stop)+1)})
        genes[gene] = {'size' : 0,
                       'low' : 0,
                       'cov' : 0}
    for line in pileup_fp:
        chrm, pos, _, cov, _, _ = line.split("\t")
        try:
            genes[chrms[chrm][pos]]['cov'] += int(cov)
            genes[chrms[chrm][pos]]['size'] += 1
            if int(cov) < 30: genes[chrms[chrm][pos]]['low'] += 1
        except KeyError:
            pass

    report_fp = open("report.txt", "w")
    report_fp.write("GENE\tSIZE\tAVG.Coverage\tLow.Positions\n")
    for gene in genes:
        report_fp.write(f"{gene[:-1]}\t{genes[gene]['size']}\t{genes[gene]['cov'] / genes[gene]['size']}\t{genes[gene]['low']}\n")
    report_fp.close()
    return 0
    
        
        
        

if __name__ == "__main__":
    bedfile_name = "/commons/Themas/Thema05/programming_challenge/CAR_0394321__en-20_target_v2.BED.txt"
    pileup_name = "/commons/Themas/Thema05/programming_challenge/pileup.txt"
    simple_pure_main(bedfile_name, pileup_name)
        
