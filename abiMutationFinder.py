#in terminal
#pip install biopython

from Bio import SeqIO

import pylab, sys, os

from Bio.Blast.Applications import NcbiblastxCommandline
        
from Bio.Blast import NCBIXML

from Bio.Alphabet import generic_dna

from Bio.Seq import Seq

import herbicideMutationFinder

path = '/home/scott/Dropbox/Programming/genomicsPython/ALS'

os.chdir('/home/scott/Dropbox/Programming/genomicsPython/ALS')


class abiAnalysis(object):

    def __init__(self, abiFile):
        self.abiFile = abiFile
        self.db = '/home/scott/Dropbox/Programming/genomicsPython/protSeqs.fasta'

    def seqObject(self):
        self.seq = SeqIO.read(self.abiFile, 'abi')

    def chromatogramData(self):
        dataG = []
        dataC = []
        dataA = []
        dataT = []
        for i in self.seq.annotations['abif_raw']['PLOC2']:
            dataG.append(self.seq.annotations['abif_raw']['DATA9'][i])
            dataA.append(self.seq.annotations['abif_raw']['DATA10'][i])
            dataT.append(self.seq.annotations['abif_raw']['DATA11'][i])
            dataC.append(self.seq.annotations['abif_raw']['DATA12'][i])

        self.dataG = dataG
        self.dataC = dataC
        self.dataA = dataA
        self.dataT = dataT

    def outputSeqs(self):
        primarySeq = ''
        secondarySeq = ''

        for i in range(len(self.dataG)):
            data = [self.dataG[i], self.dataA[i], self.dataT[i], self.dataC[i]]
            maxVal = max(data)
            minVal = min(data)
            numberMin = 0
            numberMax = 0
            primaryNuc = ''
            secondaryNuc = ''
            for value in data:
                if value < (maxVal - (0.8 * maxVal)):
                    numberMin += 1
                if value > (maxVal - (0.8 * maxVal)):
                    numberMax += 1
                if value == maxVal:
                    if data.index(value) == 0:
                        primaryNuc += 'G'
                    elif data.index(value) == 1:
                        primaryNuc += 'A'
                    elif data.index(value) == 2:
                        primaryNuc += 'T'
                    elif data.index(value) == 3:
                        primaryNuc += 'C'
                elif value > (maxVal - (0.8 * maxVal)):
                    if data.index(value) == 0:
                        secondaryNuc += 'G'
                    elif data.index(value) == 1:
                        secondaryNuc += 'A'
                    elif data.index(value) == 2:
                        secondaryNuc += 'T'
                    elif data.index(value) == 3:
                        secondaryNuc += 'C'

            if numberMax >= 3 or numberMin == 3:
                primarySeq += primaryNuc
                secondarySeq += primaryNuc
            else:
                primarySeq += primaryNuc
                secondarySeq += secondaryNuc

        self.primarySeq = primarySeq
        self.secondarySeq = secondarySeq
        
    def databaseSearchUsingLocalBlast(self):
        '''This search uses a database created using makeblastdb.  We have to create
        and store these files so they are premade for use.  
        DATABASE: PATH to databasefiles using makeblastdb.  DATABASE should include
        the base name of the files.
        QUERY: Should be the query sequence in FASTA format.  There maybe a way to 
        change this so that it is not FASTA it is just a query string if not, the 
        function will have to take the query sequence and create a text file formatted
        as FASTA before doing the search.
        OUTPUT is the user defined xml file that will be parsed later in the program.
        '''
        name = str(self.abiFile) + ".fasta"
        seq = open(name, 'w')
        seq.write(">" + str(self.abiFile))
        seq.write("\n")
        seq.write(self.primarySeq)
        seq.close()
        #seqObject = Seq(self.primarySeq, generic_dna)

        #print(self.primarySeq)
        #print(seqObject)
        #print(type(seqObject))
        
        output = "testResults.txt"
        q = name
        cline_blast = NcbiblastxCommandline(query=q, \
        db = self.db, \
        evalue = 0.00001, outfmt=5, \
        out= output)
        #Can allow for user defined evalue later.
        
        stdout, stderr = cline_blast()
        #This is the line I cannot get to run on my computer.  This line should
        #run the blast search and create the output file.  But because I am on a windows
        #machine it cannot run the commandline blast. But when I run the generated
        #commandline script it runs on ubuntu command line.  I assume this will work
        #if set up correctly.  STDERR and STDOUT will be empty.
        
        results_handle = open(output, 'r')
    
        blast_records = NCBIXML.read(results_handle)
        
        results_handle.close()
        
        for alignment in blast_records.alignments:
            print(alignment)
            for hsp in alignment.hsps:
                print(hsp)
                print(alignment.title)
                moa = alignment.title.split(" ")[0]

        results_handle.close()
        
        self.moa = moa
        
        #return blast_records.alignments

    def outputChromatogram(self):
        #pylab.figure(1)
        pylab.subplot(2, 1, 1)
        pylab.plot(self.dataG)
        pylab.plot(self.dataA)
        pylab.plot(self.dataT)
        pylab.plot(self.dataC)
        pylab.legend(("G","A","T","C"))
        pylab.subplot(2, 1, 2)
        pylab.plot(self.seq.letter_annotations['phred_quality'])
        pylab.show()


def findAB1Files(path):
    files2Analyze = []
    listDir = os.listdir(path)
    for item in listDir:
        if item.endswith("ab1"):
            files2Analyze.append(item)
    return files2Analyze
    
def batchAnalysis(path):
    ab1Files = findAB1Files(path)
    print(ab1Files)
    for ab1 in ab1Files:
        seq = abiAnalysis(ab1)
        seq.seqObject()
        seq.chromatogramData()
        seq.outputSeqs()
        seq.databaseSearchUsingLocalBlast()
        if (seq.primarySeq == seq.primarySeq) == True:
            
            stdout = sys.stdout
            nameFile = "ANALYSIS_" + str(seq.moa) + "_" + str(seq.abiFile) + "_RESULTS"
            openFile = open(nameFile, 'w')
            sys.stdout = openFile
            print("***********************************")
            print('\n')
            print(seq.abiFile)
            print('\n')
            print('***********************************')
            print('\n')
            print("Mode of Action Evaluated: " + str(seq.moa))
            print('\n')
            print('***********************************')
            print('\n')
            print("")
            print("***********************************")
            print('\n')
            print("Checking Primary Sequence")
            print('\n')
            print("***********************************")
            print('\n')
            print("")
            print('**********Primary Sequence**********')
            print('\n')
            print(seq.primarySeq)
            print('\n')
            print("************************************")
            print('\n')
            herbicideMutationFinder.mutationFinder(seq.primarySeq, seq.moa)
            sys.stdout = stdout

            openFile.close()


        else:

            stdout = sys.stdout
            nameFile = "ANALYSIS_" + str(seq.moa) + "_" + str(seq.abiFile) + "_RESULTS"
            openFile = open(nameFile, 'w')
            sys.stdout = openFile
            print("***********************************")
            print('\n')
            print(seq.abiFile)
            print('\n')
            print('***********************************')
            print('\n')
            print("Mode of Action Evaluated: " + str(seq.moa))
            print('\n')
            print('***********************************')
            print('\n')
            print("")
            print("***********************************")
            print('\n')
            print("Checking Primary Sequence")
            print('\n')
            print("***********************************")
            print('\n')
            print("")
            print('**********Primary Sequence**********')
            print('\n')
            print(seq.primarySeq)
            print('\n')
            print("************************************")
            print('\n')
            herbicideMutationFinder.mutationFinder(seq.primarySeq, seq.moa)
            print('\n')
            print("")
            print("***********************************")
            print('\n')
            print("Checking Secondary Sequence")
            print('\n')
            print("***********************************")
            print('\n')
            print("")
            print('**********Secondary Sequence**********')
            print('\n')
            print(seq.secondarySeq)
            print('\n')
            print("************************************")
            print('\n')
            herbicideMutationFinder.mutationFinder(seq.secondarySeq, seq.moa)
            print('\n')
            print("")
            sys.stdout = stdout

            openFile.close()
            
batchAnalysis(path)

#seq2 = abiAnalysis('CG34-3KU4F_KU-Acc4F_R29782_7.ab1')
#seq2 = abiAnalysis('ACCaseFastaTestCases/CG34-3KU4F_KU-Acc4F_R29782_7.ab1')
# seq2 = abiAnalysis('ALS/31-5_ALS5F_Z82221_13.ab1')
# seq2.seqObject()
# seq2.chromatogramData()
# seq2.outputSeqs()
# #seq2.outputChromatogram()
# print(seq2.primarySeq)
# print(type(seq2.primarySeq))
# seq2.databaseSearchUsingLocalBlast()
# 
# print(seq2.moa)
# stdout = sys.stdout
# nameFile = "ANALYSIS_" + str(seq2.moa) + "_" + "31-5_ALS5F_Z82221_13.ab1" + "_RESULTS"
# openFile = open(nameFile, 'w')
# sys.stdout = openFile
# print("")
# print("***********************************")
# print('\n')
# print("Checking Primary Sequence")
# print('\n')
# print("***********************************")
# print('\n')
# print("")
# 
# herbicideMutationFinder.mutationFinder(seq2.primarySeq, seq2.moa)
# 
# print(seq2.secondarySeq)
# print("")
# print("***********************************")
# print("")
# print("Checking Secondary Sequence")
# print("")
# print("***********************************")
# print("")
# herbicideMutationFinder.mutationFinder(seq2.secondarySeq, "ALS")
# 
# sys.stdout = stdout
# 
# openFile.close()

#Reading material
#https://biopython.org/wiki/SeqIO
#https://biopython.org/wiki/SeqIO#:~:text=File%20Formats%20%20%20%20Format%20name%20,determine%20th%20...%20%2028%20more%20rows%20
#https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py
#develop a function that grabs all .ab1 files
#http://www6.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf