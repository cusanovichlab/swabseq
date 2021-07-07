#!/usr/bin/env python3

import sys

def main():
    Plate_ID = 0
    #index2 = 10
    #Twist_RNA_copies = 1
    #ATCC_DNA_copies = 3
    #ATCC_virus_copies = 4
    #spike_copies = 2
    #lysate = 6
    #nCoV_amplicon = 4
    #nCoV_primer_nM = 8
    #RPP30_primer_nM = 9
    #RPP30_inner_primer_nM = 10
    #bc_set = 6
    #RT_temp = 12
    #PCR_cycles = 13
    #Sample_Well = 7
    #index = 9
    #Sample_ID = 8


    #plate_idx = 0
    #index2_idx = 10
    #index_idx = 9
    #well_idx = 7
    #rand_idx1 = -10
    #rand_idx2 = -11    

    for line in sys.stdin:
        if "Plate" not in line:
            continue
        line = line.strip()
        liner = line.split()
        data = line.split(",")
        if data[0] == "Plate_ID":
            index = data.index("index")
            index2 = data.index("index2")
            Sample_Well = data.index("Sample_Well")
            baseline = [index,index2,Sample_Well,Plate_ID]
            alls = range(len(data))
            remainder = []
            for number in alls:
                if number not in baseline:
                    remainder.append(number)
            #continue
        
        #sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(data[index] + data[index2], data[index], data[index2], data[Plate_ID], data[Sample_Well], data[nCoV_amplicon], data[lysate], data[Twist_RNA_copies], data[ATCC_RNA_copies], data[ATCC_virus_copies]))
        #sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(data[index] + data[index2], data[index], data[index2], data[Plate_ID], data[Sample_Well], data[nCoV_amplicon], data[Twist_RNA_copies], data[ATCC_DNA_copies]))
        outter = data[index] + data[index2], data[index], data[index2], data[Plate_ID], data[Sample_Well]
        outter = list(outter)
        for number in remainder:
            outter.append(data[number])
        sys.stdout.write('\t'.join(outter) + '\n')
        #sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(data[index] + data[index2], data[index], data[index2], data[Plate_ID], data[Sample_Well], data[remainder]))

    return 1

if __name__=="__main__":
    main()
