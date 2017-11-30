
'''Searches sequences to determine if specific herbicide mutations are present.
'''
import getopt, sys, os, re
from Bio.Seq import Seq
from Bio import pairwise2, SeqIO

# Get input filename from command line arguments
seqFile = sys.argv[1]
#target = sys.argv[2]
try:
    inFile = SeqIO.read(seqFile, 'abi')
except IOError:
    print('Error opening input file.')
    print('Spaces in directories can cause this error.  Use underscores.  Also try double quotes if not working.')
    sys.exit(1)

inputTarget = raw_input("For target site enter 'a' for alpha-tubulin, \n \
'b' for acetolactate synthase, \n \
'c' for acetylcoa carboxylase, \n \
'd' for protoporphyrinogen oxidase, \n \
'e' for epsp synthase, \n \
'f' for psbA PSII D1 protein. \n \
Enter a selection: ")

if inputTarget == 'a':
    target = 'alpha-tubulin'
elif inputTarget == 'b':
    target = 'ALS'
elif inputTarget == 'c':
    target = 'ACCase'
elif inputTarget == 'd':
    target = 'PPO'
elif inputTarget == 'e':
    target = 'EPSP'
elif inputTarget == 'f':
    target = 'psbA'
else:
    print "Error entering target site, must enter a, b, c, d or e"
    sys.exit(1)
    
# seqFile = ''
# target = ''
# species_input = ''
# try:
#     opts, args = getopt.getopt(sys.argv, "hf:t:", ["ab1=", "target="])
# except getopt.GetoptError:
#     print('herbicideMutationFinder<VERSION>.py -f <inputAB1Files> -t <targetSite>')
#     sys.exit(2)
# for opt, arg in opts:
#     if opt == '-h':
#         print('herbicideMutationFinder<VERSION>.py -f <inputAB1Files> -t <targetSite>')
#         sys.exit()
#     elif opt in ("-f", "--ab1"):
#         seqFile = arg
#     elif opt in ("-t", "--target"):
#         target = arg


# Iterate through file, add items to database
# Open file
# try:
#     infile = SeqIO.read(seqFile, 'abi')
# except IOError:
#     print('Error opening input file.')
#     print('Spaces in directories can cause this error.  Use underscores.  Also try double quotes if not working.')
#     sys.exit(1)

class resistanceMutationSearch(object):
    
    def __init__(self, sequence, targetSite):
        '''Initiate the object with a sequence and a target-site.  Target-site 
        can be selected by user from dropdown menu or provided by raw input.
        '''
        self.sequence = sequence
        
        self.targetSite = targetSite
           
        
    def determineTarget(self):
        '''
        DETERMINETARGET:  Based on input targetSite, sets self.segments, 
        self.susSegments (motifs surrounding susceptible target sites), 
        self.resSegments (a dict of susceptible and resistant motifs as key value pairs),
        self.compDNA (string of comparision DNA for alignments)
        self.compProt (string of comparison Prot for alignments)
        These items will be used in subsequent functions to find mutations.
        '''
        if self.targetSite == 'ALS':
            self.segments = ['MVVQ', 'ELDQQ', 'PENES', 'GGASM', 'QVPR', 'TDAF','FDDR','QWED', 'PSGG', 'GLPL']
            self.susSegments = ['GGASM', 'QVPR', 'TDAF','FDDR','QWED', 'PSGG']
            self.resSegments = {'GGASM':'GG\wSM','QVPR':'QV\wR','TDAF':'TD\wF', 'FDDR':'FD\wR',"QWED":"Q\wED", 'PSGG':'P\w\wG'}
            self.compDNA = 'GCCACAGCCACAGCCACGTCCACAGCCGTCGCCATCTCGGGCGCCACCTCCGCCCTACCCAAACCCACCC\
            TCCCGCGCCACCTGCCCGCCCCACGCCGCGCCGCCCTCGCCGCCGCCACCCGCATCAGGTGCTCCACGGT\
            GTCCCCTTCGCCCGCCCCTCCCGCCACCGCGCTCCGCCCATGGGGCCCCACCGAGCCCCGCAAGGGCGCC\
            GACATCCTCGTCGAGGCCCTGGAGCGCTGCGGCATCAGCGACGTCTTCGCCTACCCCGGCGGCGCCTCGA\
            TGGAGATCCACCAGGCGCTCACGCGCTCGCCGGCCATCACCAACCACCTCTTCCGGCACGAGCAGGGGGA\
            GGCGTTCGCCGCGTCCGGGTACGCCCGCGCCTCCGGCCGCGTCGGGGTCTGCGTCGCCACCTCCGGCCCC\
            GGCGCCACCAACCTCGTCTCCGCGCTCGCCGACGCTCTGCTCGACTCCATCCCGATGGTCGCCATCACGG\
            GGCAGGTCCCGCGCCGCATGATCGGCACGGACGCCTTCCAGGAGACGCCGATTGTGGAGGTCACCCGTTC\
            CATCACCAAGCACAATTACCTGGTCCTTGACGTGGAGGACATCCCCCGCGTCATTCAGGAAGCCTTCTTC\
            CTCGCCTCCTCCGGCCGGCCGGGGCCGGTGCTGGTCGACATCCCCAAGGACATCCAGCAGCAGATGGCCG\
            TGCCTGTCTGGGACGCGCCAATGAGTCTGCCAGGGTACATTGCTCGCCTCCCTAAGCCGCCGGCTACCGA\
            ATTGCTTGAGCAGGTCCTGCGTCTGGTTGGTGAGGCTCGGCGCCCAATTCTGTATGTTGGTGGTGGCTGC\
            TCTGCGTCCGGCGAGGAGTTGCGCCGCTTTGTTGAGCTCACTGGGATCCCAGTGACAACTACCCTCATGG\
            GTCTTGGCAACTTCCCCAGCGATGACCCACTGTCTCTGCGTATGCTTGGGATGCATGGTACAGTCTACGC\
            CAATTACGCGGTAGATAAGGCTGACCTGCTGCTTGCATTTGGTGTGCGGTTTGATGACCGTGTGACTGGA\
            AAAATAGAGGCTTTTGCAAGCAGGTCCAAGATTGTGCACATTGACATTGATCCAGCTGAGATTGGCAAGA\
            ACAAGCAGCCACACGTCTCCATTTGTGCAGATGTCAAGATCGCTTTGGAGGGCTTGAATTCTCTTCTGCT\
            AAATGGGAGCAAAACACACAAGAGTTTAGATTTTAGTTCGTGGCATGAGGAGTTGGACCAGCAGAAGAGG\
            GAGTTTCCTCTGGGATTCAAAACTTTTGGTGAGGCGATCCCACCACAATATGCTATCCAGGTACTGGATG\
            AGCTGACCAAAGGGGAGGCGATCATTGCCACTGGTGTTGGGCAGCACCAGATGTGGGCGGCTCAGTATTA\
            CACGTACAAGCGGCCACGTCAGTGGCTGTCTTCGGCTGGTCTTGGAGCAATGGGGTTTGGGTTGCCAGCT\
            GCAGCTGGTGCTGCTGTGGCCAACCCAGGTGTCACAGTTGTTGACATTGATGGAGATGGTAGCTTCCTCA\
            TGAATATTCAGGAGTTGGCACTGATTCGTATTGAGAACCTCCCTGTTAAGGTGATGATACTGAACAACCA\
            ACATCTGGGAATGGTGGTGCAGTGGGAGGACAGGTTTTACAAGGCCAATCGGGCGCACACGTACCTTGGG\
            AACCCAGAAAATGAGAGTGAGATATATCCAGATTTTGTGACGATTGCCAAGGGGTTCAATGTTCCTGCTG\
            TTCGTGTGACAAAGAAAAGTGAAGTCCGTGCAGCAATCAAGAAGATGCTTGAGACTCCAGGGCCATACTT\
            GTTGGATATCATCGTCCCTCACCAGGAGCATGTGCTGCCTATGATCCCCAGCGGTGGTGCTTTCAAGGAC\
            ATTATCATGGAGGGCGATGGCA'
            self.compProt = 'MATATATSTAVAISGATSALPKPTLPRHLPAPRRAALAAATRIRC\
            STVSPSPAPPATALRPWGPTEPRKGADILVEALERCGISDVFAYPGGASMEIHQALTR\
            SPAITNHLFRHEQGEAFAASGYARASGRVGVCVATSGPGATNLVSALADALLDSIPMV\
            AITGQVPRRMIGTDAFQETPIVEVTRSITKHNYLVLDVEDIPRVIQEAFFLASSGRPG\
            PVLVDIPKDIQQQMAVPVWDAPMSLPGYIARLPKPPATELLEQVLRLVGEARRPILYV\
            GGGCSASGEELRRFVELTGIPVTTTLMGLGNFPSDDPLSLRMLGMHGTVYANYAVDKA\
            DLLLAFGVRFDDRVTGKIEAFASRSKIVHIDIDPAEIGKNKQPHVSICADVKIALEGL\
            NSLLLNGSKTHKSLDFSSWHEELDQQKREFPLGFKTFGEAIPPQYAIQVLDELTKGEA\
            IIATGVGQHQMWAAQYYTYKRPRQWLSSAGLGAMGFGLPAAAGAAVANPGVTVVDIDG\
            DGSFLMNIQELALIRIENLPVKVMILNNQHLGMVVQWEDRFYKANRAHTYLGNPENES\
            EIYPDFVTIAKGFNVPAVRVTKKSEVRAAIKKMLETPGPYLLDIIVPHQEHVLPMIPSGGAFKDIIMEGD'
            self.testAccession = 'KM388812'
        if self.targetSite == "ACCase":
            self.segments = ["GREV", "PDDL", "IPENT", "VENIH", "GQVWFP", "LANWRG", "LANW","FEGIL", "LFEGILQ", "IDSKIN", "TAKGN"]
            self.susSegments = ["VENIH", "GQVWFP", "LANWRG", "LFEGILQ", "IDSKIN"]
            self.resSegments = {"VENIH":"VEN\wH", "GQVWFP":"GQV\wFP", "LANWRG":"LAN\wRG", "LFEGILQ":"LFEG\wLQ", "IDSKIN":"I\wSKIN", "TAKGN":"TAK\wN"}
            self.compDNA = 'GCCACAGCCACAGCCACGTCCACAGCCGTCGCCATCTCGGGCGCCACCTCCGCCCTACCCAAACCCACCC\
            TCCCGCGCCACCTGCCCGCCCCACGCCGCGCCGCCCTCGCCGCCGCCACCCGCATCAGGTGCTCCACGGT\
            GTCCCCTTCGCCCGCCCCTCCCGCCACCGCGCTCCGCCCATGGGGCCCCACCGAGCCCCGCAAGGGCGCC\
            GACATCCTCGTCGAGGCCCTGGAGCGCTGCGGCATCAGCGACGTCTTCGCCTACCCCGGCGGCGCCTCGA\
            TGGAGATCCACCAGGCGCTCACGCGCTCGCCGGCCATCACCAACCACCTCTTCCGGCACGAGCAGGGGGA\
            GGCGTTCGCCGCGTCCGGGTACGCCCGCGCCTCCGGCCGCGTCGGGGTCTGCGTCGCCACCTCCGGCCCC\
            GGCGCCACCAACCTCGTCTCCGCGCTCGCCGACGCTCTGCTCGACTCCATCCCGATGGTCGCCATCACGG\
            GGCAGGTCCCGCGCCGCATGATCGGCACGGACGCCTTCCAGGAGACGCCGATTGTGGAGGTCACCCGTTC\
            CATCACCAAGCACAATTACCTGGTCCTTGACGTGGAGGACATCCCCCGCGTCATTCAGGAAGCCTTCTTC\
            CTCGCCTCCTCCGGCCGGCCGGGGCCGGTGCTGGTCGACATCCCCAAGGACATCCAGCAGCAGATGGCCG\
            TGCCTGTCTGGGACGCGCCAATGAGTCTGCCAGGGTACATTGCTCGCCTCCCTAAGCCGCCGGCTACCGA\
            ATTGCTTGAGCAGGTCCTGCGTCTGGTTGGTGAGGCTCGGCGCCCAATTCTGTATGTTGGTGGTGGCTGC\
            TCTGCGTCCGGCGAGGAGTTGCGCCGCTTTGTTGAGCTCACTGGGATCCCAGTGACAACTACCCTCATGG\
            GTCTTGGCAACTTCCCCAGCGATGACCCACTGTCTCTGCGTATGCTTGGGATGCATGGTACAGTCTACGC\
            CAATTACGCGGTAGATAAGGCTGACCTGCTGCTTGCATTTGGTGTGCGGTTTGATGACCGTGTGACTGGA\
            AAAATAGAGGCTTTTGCAAGCAGGTCCAAGATTGTGCACATTGACATTGATCCAGCTGAGATTGGCAAGA\
            ACAAGCAGCCACACGTCTCCATTTGTGCAGATGTCAAGATCGCTTTGGAGGGCTTGAATTCTCTTCTGCT\
            AAATGGGAGCAAAACACACAAGAGTTTAGATTTTAGTTCGTGGCATGAGGAGTTGGACCAGCAGAAGAGG\
            GAGTTTCCTCTGGGATTCAAAACTTTTGGTGAGGCGATCCCACCACAATATGCTATCCAGGTACTGGATG\
            AGCTGACCAAAGGGGAGGCGATCATTGCCACTGGTGTTGGGCAGCACCAGATGTGGGCGGCTCAGTATTA\
            CACGTACAAGCGGCCACGTCAGTGGCTGTCTTCGGCTGGTCTTGGAGCAATGGGGTTTGGGTTGCCAGCT\
            GCAGCTGGTGCTGCTGTGGCCAACCCAGGTGTCACAGTTGTTGACATTGATGGAGATGGTAGCTTCCTCA\
            TGAATATTCAGGAGTTGGCACTGATTCGTATTGAGAACCTCCCTGTTAAGGTGATGATACTGAACAACCA\
            ACATCTGGGAATGGTGGTGCAGTGGGAGGACAGGTTTTACAAGGCCAATCGGGCGCACACGTACCTTGGG\
            AACCCAGAAAATGAGAGTGAGATATATCCAGATTTTGTGACGATTGCCAAGGGGTTCAATGTTCCTGCTG\
            TTCGTGTGACAAAGAAAAGTGAAGTCCGTGCAGCAATCAAGAAGATGCTTGAGACTCCAGGGCCATACTT\
            GTTGGATATCATCGTCCCTCACCAGGAGCATGTGCTGCCTATGATCCCCAGCGGTGGTGCTTTCAAGGAC\
            ATTATCATGGAGGGCGATGGCA'
            self.compProt = 'NPERGFKYIYLNEEDYGRISSSVIAHKTQLDSGEIRWVIDSVVG\
            KEDGLGVENIHGSAAIASAYSRAYEETFTLTFVSGRTVGIGAYLARLGIRCIQRIDQP\
            IILTGFSALNKLLGREVYSSHMQLGGPKIMATNGVVHLTVPDDLEGVSNILRWLSYVP\
            ANIGGPLPITKSLDPIDRPVAYIPENTCDPRAAISGIXDXQGKWLGGMFDKDSFVETF\
            EGWAKTVVTGRAKLGGIPVGVIAVETQTMMQLVPADPGQPDSHERSVPRAGQVWFPDS\
            ATKTAQAMLDFNREGLPLFILANWRGFSGGQRDLFEGILQAGSTIVENLRTYNQPAFV\
            YIPKAAELRGGAWVVIDSKINPDRIECYAETTAKGNVLEPQGLIEIKFRSEXLQDCMD\
            RLDPELVNLKAKLQGAKHEN'
            self.testAccession = 'KJ606972'
        if self.targetSite == "psbA":
            self.segments = ['GSLV', 'QYAS','RETTENES']
            self.susSegments = ['GSLV', 'QYAS']
            self.resSegments = {'GSLV': 'GSL\w', 'QYAS':'QYA\w'} # V and S are mutation sites
            self.compDNA = 'ATGATCCCTACCTTATTGACTGCAACTTCTGTATTTATTATAGCCTTCATAGCTGCTCCTCCAGTAGATA\
            TTGATGGTATTCGTGAACCTGTTTCTGGATCTCTACTTTACGGAAACAATATTATTTCGGGTGCTATTAT\
            TCCTACTTCTGCAGCTATTGGGTTGCACTTTTACCCAATCTGGGAAGCGGCATCAGTTGATGAGTGGTTA\
            TACAATGGTGGTCCTTATGAACTAATCGTTCTACACTTCTTACTTGGTGTAGCTTGTTATATGGGTCGTG\
            AGTGGGAACTTAGTTTCCGTCTGGGTATGCGTCCGTGGATTGCTGTTGCATATTCAGCTCCGGTTGCAGC\
            GGCTACTGCTGTTTTCTTGATCTACCCAATCGGTCAAGGAAGCTTTTCTGATGGTATGCCTCTAGGAATC\
            TCTGGTACTTTCAACTTTATGATCGTATTCCAGGCTGAGCACAACATCCTTATGCACCCATTTCACATGT\
            TAGGTGTAGCTGGTGTATTCGGCGGCTCCCTATTTAGTGCTATGCATGGTTCCTTGGTAACTTCTAGTTT\
            GATCAGGGAAACCACAGAAAATGAATCTGCTAACGAAGGTTACAGATTCGGTCAAGAGGAAGAAACTTAT\
            AACATCGTAGCTGCTCATGGTTATTTTGGTCGATTGATCTTCCAATATGCTAGTTTCAACAACTCTCGTT\
            CTTTACACTTCTTCTTAGCTGCTTGGCCGGTAATCGGTATTTGGTTTACTGCTTTGGGTATTAGTACTAT\
            GGCTTTCAACCTAAACGGTTTCAACTTCAACCAATCTGTAGTTGATAGTCAAGGTCGTGTAATTAACACC\
            TGGGCTGATATCATTAACCGTGCTAACCTTGGTATGGAAGTTATGCATGAACGTAATGCTCATAACTTCC\
            CTCTAGACTTAGCTGCTATCGAAGCTCCATCTACAAATGGATAA'
            self.compProt = 'MIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISG\
            AIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLG\
            MRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHP\
            FHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYRFGQEEETYNIVAAHGY\
            FGRLIFQYASFNNSRSLHFFLAAWPVIGIWFTALGISTMAFNLNGFNFNQSVVDSQGR\
            VINTWADIINRANLGMEVMHERNAHNFPLDLAAIEAPSTNG'
            self.testAccession = 'K01200'
        if self.targetSite == "EPSPS":
            self.segments = ['KEEVQ', 'TAAGG', 'MRPL', 'LGNAG' ]
            self.susSegments = ['MRPL', 'LGNAG']
            self.resSegments = {'MRPL':'MR\wL', 'LGNAG':'LGNA\w'}
            self.compDNA = 'GCGGAGGAGGTGGTGCTGCAGCCCATCAAGGAGATCTCCGGCGTCGTGAAGCTGCCGGGGTCCAAGTCGC\
            TCTCCAACCGGATCCTCCTGCTCTCCGCCCTCGCCGAGGTAAGAAGAAGGATCCCCCCTCCCTTTCAGAG\
            TCAATTGAAATTGGATGTGGAGATGAGATTTTACCAGGGGTTAGGTGATGATTTCCTGCTGCTGAAATGT\
            CTGTGTAGTCATAGGAATATACGGATCAGATGGAGCTCGGTCATTCACTAGTAGGCATTCGTCTCTGTAG\
            ATTGTATTCACCTTAATTAAGGCGATTGATGAGGTAGACTGAAAGCAGCAGGGGGATTTGGCTATTTAGC\
            GATACTGAATCGTTGGGATGGCATGTGTTCTGTATCAGGGTACTAATTGCAATACGTGCTGAAGGTTCCA\
            TCCAGGTTGGAACTTGATCACGAGAGCCCACATTTTAATTGAAAAATGTTGACTGGATTATAGTTTCATG\
            TGTCACTGTTTATTTGGCTTCCATTAGTTGGGAAAAACTAGTCGAGGTGGAGCGTGCATTTTTAAACCCT\
            TATACCAAGATAAGAAGCAGCCGCCTTCTCCTTTTCGTTTCCTCTTAGAAAGGGGGAATTGGAACATTAG\
            ATGAGTAACGGTTGATGCTGTAACCTTTTCTCAGGGAACAACTGTGGTGGATAACCTTTTAAACAGTGAG\
            GACGTCCACTACATGCTCGGGGCCCTGAAAACCCTCGGACTCTCTGTGGAAGCGGACAAAGCTGCCAAAA\
            GAGCGGTAGTTGTTGGCTGTGGTGGCAAGTTCCCAGTTGAGAAGGATGCGAAAGAGGAGGTGCAGCTCTT\
            CTTGGGGAATGCTGGAACTGCAATGCGACCATTGACAGCAGCCGTAACTGCTGCTGGAGGAAATGCAACG\
            TGAGTTGGTTTTTCCATCCTCAGAATATGCCCGTGGAACTGAGTAGCGAAATTGTGGTGATATTTCGTGA\
            CTTATCGTGCATCTTTTCTGAATTCCAGTTATGTGCTTGATGGAGTGCCAAGAATGCGGGAGAGACCCAT\
            TGGCGACTTGGTTGTCGGATTGAAACAGCTTGGTGCAGATGTTGATTGTTTCCTTGGCACTGACTGCCCA\
            CCTGTTCGTGTCAAGGGAATCGGAGGGCTACCTGGTGGCAAGGTTAGTTACTAAAACCAGCATGCTACAT\
            TCTTCTGTAACCAATTGAAATTTTCTTGTTGAGCTTGTGCATTTCAAAGAAAAAGTCTAGAATATGTGTT\
            TTGACTTTATAAACCGCACAGATTGTCATCTCTGACCAGGAATCCATATTTTCAATAAACTCATCAAGTA\
            CTATTTGGCTCATCATCTTGCAAGTTATTTGCTCCTTTTGGGCTGGTGTTAGTTTTCTGGTCTCACATTG\
            TGTGGGATACCATGAATTAGCTCTTAAAAATGTGTTACTTTCTACAGGTTAAGTTATCTGGTTCTATCAG\
            CAGTCAGTACTTGAGTGCCTTGCTGATGGCTGCTCCTTTAGCTCTTGGGGATGTGGAGATTGAAATMATT\
            GATAAACTGATCTCCATCCCTTATGTTGAAATGACATTGAGATTGATGGAGCGTTTTGGCGTGAAAGCAG\
            AGCATTCTGATAGCTGGGACAGATTCTACATCAAGGGAGGTCAAAAATACAAGTAAGCCAGTCATTTTGT\
            TCTCAGCTTACCTTCCTAATCACTGATTTATCCCAGATAAACTGGAGACTGAAATTATTTCAGAAAAAAA\
            ACTTAGTATTGATTATGTAAGGTATTTAATCAGCTAGAATATGACCAGTCTGTTGGAGTCATGAAAGCTA\
            ACTAAAATAGCTAAGCAATGTCTAGTAAGTCAATCCATTTGGTTATAAATAGAGTTGATAACACCAAAAG\
            TAGCAACCATTCTAGGTGAGGACTTCCTGCTGCATGGCATGGGTAGCATTTTTTTTAGGGAAAGGTACCA\
            GATTAGTTTCTTTGCACCCCCTTGATGTGAACCCGCWCATGAAGTATCAGGATAACTGTCCTTATTTTGA\
            CTTAAAATACTCTCTTAGTCTTGTGTACTCAGATTGGTAGTTCAAACTTACCAGGTCCCCTAAAAATGCC\
            TACGTGGAAGGTGATGCCTCAAGTGCGAGCTATTTCTTGGCTGGTGCTGCAATCACTGGAGGGACTGTGA\
            CTGTTGAAGGTTGTGGCACCACCAGTCTGCAGGTAGAAATCTGCCAAAATGCTACTTTTGATTCATCTGC\
            CATAATACTCACAATAGACATATTTGATTCACCTGCCAATGCCAAGTGATGTGTTAATGAGTGGATACTC\
            ATGTAAAGGGGCTAAAATTCTGAAACCAACAATTCTATAAATGACTAAATAAATATGCTTTCTACTTGGT\
            CGTATTCGTTGACAAATATGCCTAATGTTGCAGGGTGATGTGAAATTTGCCGAGGTACTCGAGATGATGG\
            GAGCGAAGGTTACATGGACTGAAACTAGCGTAACTGTTACCGGTCCACAACGTGAGCCATTTGGGAGGAA\
            ACACCTAAAAGCTATTGATGTTAACATGAACAAAATGCCCGATGTCGCCATGACTCTTGCCGTGGTTGCC\
            CTATTTGCTGATGGCCCAACTGCTATCAGAGATGGTAAACATCCTCAGCCGAACAATGCCTGACTTCATA\
            GCACTTTAAGCTGGTTTTGGCTATGCTTACTGAAACTTCTGTTCTTTAACAGTGGCTTCCTGGAGAGTAA\
            AGGAGACCGAGAGGATGGTTGCAATCCGGACTGAGCTAACAAAGGTAAGAATACTTCATAAACCTTTCGT\
            TTATTGCCAGATTTGGACATCTCTTTTGCTCTCCTGCCGGCGTTCCATCTTCAGCTGTTCACATGCTGAC\
            GGCTATTTACATCTTTCTAGCTGGGAGCGTCGGTCGAGGAAGGACCGGACTACTGCATTATCACACCGCC\
            CGAGAAGCTGAACGTAACGGCCATCGACACCTACGATGACCACAGGATGGCCATGGCCTTCTCCCTCGCC\
            GCCTGCGCCGACGTGCCTGTGACCATCCGGGACCCCGGCTGCACCCGCAAGACCTTCCCAGACTACTTCG\
            ACGTGCTGAGCACTTTCGTCAAGAAC'
            self.compProt = "AEEVVLQPIKEISGVVKLPGSKSLSNRILLLSALAEGTTVVDNL\
            LNSEDVHYMLGALKTLGLSVEADKAAKRAVVVGCGGKFPVEKDAKEEVQLFLGNAGTA\
            MRSLTAAVTAAGGNATYVLDGVPRMRERPIGDLVVGLKQLGADVDCFLGTDCPPVRVK\
            GIGGLPGGKVKLSGSISSQYLSALLMAAPLALGDVEIEIIDKLISIPYVEMTLRLMER\
            FGVKAEHSDSWDRFYIKGGQKYKSPKNAYVEGDASSASYFLAGAAITGGTVTVEGCGT\
            TSLQGDVKFAEVLEMMGAKVTWTETSVTVTGPQREPFGRKHLKAIDVNMNKMPDVAMT\
            LAVVALFADGPTAIRDVASWRVKETERMVAIRTELTKLGASVEEGPDYCIITPPEKLN\
            VTAIDTYDDHRMAMAFSLAACADVPVTIRDPGCTRKTFPDYFDVLSTFVKN"
            #change MRSL to MRPL for non-resistant standard
            self.testAccession = 'JN004269'#Contains a S-106 mutation so changed to P
        if self.targetSite == "alpha-tubulin" :
            self.segments = ['NEFQT', 'SLTAS', 'HFMLS', 'SAEKA', 'DVNEF']
            self.susSegments = ['SLTAS', 'HFMLS']
            self.resSegments = {'SLTAS': 'SL\wAS', 'HFMLS': 'HF\wLS'} #Thr-239-Ile, Met-268-Thr, Leu-136-Phe
            self.compDNA = 'CAGAGACAGGCGTCTTCGTACTCGCCTCTCTCCGCGACTCCAAGCTTTCTCCCTCCTCCCATTTCCCGTC\
            GCCGCCGCCTCACCCGCCCGACACCATGAGGGAGTGCATCTCGATCCACATCGGCCAGGCCGGTATCCAG\
            GTCGGAAACGCTTGCTGGGAGCTCTACTGCCTCGAGCATGGCATCCAGGCTGACGGTCAGATGCCCGGTG\
            ACAAGACCATTGGAGGAGGTGATGATGCTTTCAACACCTTCTTCAGTGAGACTGGCGCCGGCAAGCATGT\
            GCCCCGTGCCGTGTTTGTTGACCTTGAGCCCACTGTGATCGATGAGGTCAGGACTGGCACCTACCGCCAG\
            CTGTTTCACCCTGAGCAGCTCATCAGTGGCAAGGAGGATGCTGCCAACAACTTTGCCCGTGGTCACTACA\
            CCATTGGCAAGGAGATTGTTGACCTGTGCCTTGACCGCATCAGGAAGCTTGCCGACAACTGTACTGGTCT\
            CCAGGGCTTCCTTGTCTTCAACGCTGTCGGTGGAGGAACGGGCTCTGGTCTTGGTTCCCTCCTCCTTGAG\
            CGCCTGTCTGTTGACTACGGCAAGAAGTCCAAGCTCGGGTTCACTGTCTACCCGTCTCCTCAGGTCTCCA\
            CCTCGGTGGTTGAGCCATACAACAGTGTGCTGTCCACCCACTCCCTCCTTGAGCACACCGATGTGGCTGT\
            GCTGCTTGACAACGAGGCCATCTACGACATCTGCCGCCGCTCCCTGGACATTGAGCGCCCAACCTACACC\
            AACCTGAACAGGCTTGTTTCTCAGGTCATTTCATCACTGACAGCCTCTCTGAGGTTCGATGGTGCTCTGA\
            ACGTGGATGTGAACGAGTTCCAGACCAACTTGGTGCCCTACCCGAGGATCCACTTCATGCTTTCATCCTA\
            CGCTCCAGTGATCTCCGCGGAGAAGGCCTACCACGAGCAGCTGTCCGTGGCTGAGATCACCAACAGCGCG\
            TTCGAGCCTTCCTCCATGATGGCCAAGTGCGACCCCCGCCACGGCAAGTACATGGCCTGCTGCCTCATGT\
            ACCGTGGTGATGTGGTGCCCAAGGACGTGAACGCCGCCGTTGCCACCATCAAGACCAAGCGCACCATCCA\
            GTTCGTGGACTGGTGCCCCACTGGCTTCAAGTGCGGTATCAACTACCAGCCACCCAGCGTCGTCCCCGGC\
            GGCGACCTGGCCAAGGTGCAGAGGGCCGTGTGCATGATCTCCAACTCCACCAGTGTCGTCGAGGTGTTCT\
            CCCGCATCGACCACAAGTTCGACCTCATGTACGCCAAGCGCGCCTTCGTCCACTGGTACGTGGGTGAGGG\
            TATGGAGGAGGGTGAGTTCTCCGAGGCGCGTGAGGACCTTGCTGCCCTTGAGAAGGACTACGAGGAGGTC\
            GGCGCTGAGTTCGACGAGGGTGAGGAAGGTGATGAGGGTGACGAGTACTAGATGAATCTACCGCTTCCTG\
            CTGTTGTGTCAGGCCTGTGTGCCGCTGCTATCCTGTGATCTGCCCGAGGGCGCTATCGTGTCGTGTCAGT\
            TTGAACTATTTGTCATTGTGTGGTTACAACCCCTGAAGTTGTAGACATGTTTAATTCCCCGCTTTGCTAC\
            TGGGTTATCAACATCGTTATGTTTGTCTAGATTGCTGCCCCTTTGTTTTGTTTCTTTTTTGTGCGTGGTT\
            CTTCTCTTTTATC'
            self.compProt = 'MRECISIHIGQAGIQVGNACWELYCLEHGIQADGQMPGDKTIGG\
GDDAFNTFFSETGAGKHVPRAVFVDLEPTVIDEVRTGTYRQLFHPEQLISGKEDAANN\
FARGHYTIGKEIVDLCLDRIRKLADNCTGLQGFLVFNAVGGGTGSGLGSLLLERLSVD\
YGKKSKLGFTVYPSPQVSTSVVEPYNSVLSTHSLLEHTDVAVLLDNEAIYDICRRSLD\
IERPTYTNLNRLVSQVISSLTASLRFDGALNVDVNEFQTNLVPYPRIHFMLSSYAPVI\
SAEKAYHEQLSVAEITNSAFEPSSMMAKCDPRHGKYMACCLMYRGDVVPKDVNAAVAT\
IKTKRTIQFVDWCPTGFKCGINYQPPSVVPGGDLAKVQRAVCMISNSTSVVEVFSRID\
HKFDLMYAKRAFVHWYVGEGMEEGEFSEAREDLAALEKDYEEVGAEFDEGEEGDEGDEY'
            self.testAccession = 'AF008120'
        if self.targetSite == "PPO": #G210 deletion, R98G, R98M
            self.segments = ['VAGTCGGDP', 'QNKRYI', 'VSTKN', 'KLKSHG', 'TFPEV', 'GGENAS', 'TAPIRN']
            self.susSegments = ['GTCGGDP', 'QNKRYI']
            self.resSegments = {'GTCGGDP':'GTC\w?GDP', 'QNKRYI': 'QNK\wYI'} 
            self.compDNA = 'ATGGTAATTCAATCCATTACCCACCTTTCACCAAACCTTGCATTGCCATCGCCATTGTCAGTTTCAACCA\
AGAACTACCCAGTAGCTGTAATGGGCAACATTTCTGAGCGGGAAGAACCCACTTCTGCTAAAAGGGTTGC\
TGTTGTTGGTGCTGGAGTTAGTGGACTTGCTGCTGCATATAAGCTAAAATCCCATGGTTTGAGTGTGACA\
TTGTTTGAAGCTGATTCTAGAGCTGGAGGCAAACTTAAAACTGTTAAAAAAGATGGTTTTATTTGGGATG\
AGGGGGCAAATACTATGACAGAAAGTGAGGCAGAGGTCTCGAGTTTGATCGATGATCTTGGGCTTCGTGA\
GAAGCAACAGTTGCCAATTTCACAAAATAAAAGATACATAGCTAGAGCCGGTCTTCCTGTGCTACTACCT\
TCAAATCCCGCTGCACTACTCACGAGCAATATCCTTTCAGCAAAATCAAAGCTGCAAATTATGTTGGAAC\
CATTTCTCTGGAGAAAACACAATGCTACTGAACTTTCTGATGAGCATGTTCAGGAAAGCGTTGGTGAATT\
TTTTGAGCGACATTTTGGGAAAGAGTTTGTTGATTATGTTATTGACCCTTTTGTTGCGGGTACATGTGGT\
GGAGATCCTCAATCGCTTTCCATGCACCATACATTTCCAGAAGTATGGAATATTGAAAAAAGGTTTGGCT\
CTGTGTTTGCCGGACTAATTCAATCAACATTGTTATCTAAGAAGGAAAAGGGTGGAGAAAATGCTTCTAT\
TAAGAAGCCTCGTGTACGTGGTTCATTTTCATTTCAAGGTGGAATGCAGACACTTGTTGACACAATGTGC\
AAACAGCTTGGTGAAGATGAACTCAAACTCCAGTGTGAGGTGCTGTCCTTGTCATATAACCAGAAGGGGA\
TCCCCTCACTAGGGAATTGGTCAGTCTCTTCTATGTCAAATAATACCAGTGAAGATCAATCTTATGATGC\
TGTGGTTGTCACTGCTCCAATTCGCAATGTCAAAGAAATGAAGATTATGAAATTTGGAAATCCATTTTCA\
CTTGACTTTATTCCAGAGGTGACGTACGTACCCCTTTCCGTTATGATTACTGCATTCAAAAAGGATAAAG\
TGAAGAGACCTCTTGAGGGCTTCGGAGTTCTTATCCCCTCTAAAGAGCAACATAATGGACTGAAGACTCT\
TGGTACTTTATTTTCCTCCATGATGTTTCCTGATCGTGCTCCATCTGACATGTGTCTCTTTACTACATTT\
GTCGGAGGAAGCAGAAATAGAAAACTTGCAAACGCTTCAACGGATGAATTGAAGCAAATAGTTTCTTCTG\
ACCTTCAGCAGCTGTTGGGCACTGAGGACGAACCTTCATTTGTCAATCATCTCTTTTGGAGCAACGCATT\
CCCATTGTATGGACACAATTACGATTCTGTTTTGAGAGCCATAGACAAGATGGAAAAGGATCTTCCTGGA\
TTTTTTTATGCAGGTAACCATAAGGGTGGACTTTCAGTGGGAAAAGCGATGGCCTCCGGATGCAAGGCTG\
CGGAACTTGTAATATCCTATCTGGACTCTCATATATACGTGAAGATGGATGAGAAGACCGCGTAA'
            self.compProt = 'MVIQSITHLSPNLALPSPLSVSTKNYPVAVMGNISEREEPTSAK\
RVAVVGAGVSGLAAAYKLKSHGLSVTLFEADSRAGGKLKTVKKDGFIWDEGANTMTES\
EAEVSSLIDDLGLREKQQLPISQNKRYIARAGLPVLLPSNPAALLTSNILSAKSKLQI\
MLEPFLWRKHNATELSDEHVQESVGEFFERHFGKEFVDYVIDPFVAGTCGGDPQSLSM\
HHTFPEVWNIEKRFGSVFAGLIQSTLLSKKEKGGENASIKKPRVRGSFSFQGGMQTLV\
DTMCKQLGEDELKLQCEVLSLSYNQKGIPSLGNWSVSSMSNNTSEDQSYDAVVVTAPI\
RNVKEMKIMKFGNPFSLDFIPEVTYVPLSVMITAFKKDKVKRPLEGFGVLIPSKEQHN\
GLKTLGTLFSSMMFPDRAPSDMCLFTTFVGGSRNRKLANASTDELKQIVSSDLQQLLG\
TEDEPSFVNHLFWSNAFPLYGHNYDSVLRAIDKMEKDLPGFFYAGNHKGGLSVGKAMA\
SGCKAAELVISYLDSHIYVKMDEKTA'   
            self.testAccession = 'DQ386117'                    
            
    def idCorrectProtein(self):
        possibleProts = self.translator()
        self.prot = None
        for prot in possibleProts.values():
            for seg in self.segments:
                if seg in prot:
                    self.prot = prot
                    break
    
    def alignProteins(self):
        align = pairwise2.align.globalxx(self.compProt, self.prot)
        formattedAlignment = (pairwise2.format_alignment(*align[0]))
        
    def detectMutations(self):
        '''Determine if a specific mutation is present in the the selected
        protein string.  Return DICTIONARY of mutations found and a DICTIONARY 
        of Normal Segments and Mutant Segments
        '''
        mutationsPresent = {}
        specificMutations = {}
        for seg in self.susSegments:
            if seg in self.prot:
                mutationsPresent[seg]="No Mutation Present"
            elif re.search(self.resSegments[seg], self.prot):
                reObj = re.search(self.resSegments[seg], self.prot)
                mutationsPresent[seg]="Mutation Present"
                start = reObj.start() #self.resSegments[seg].find(self.prot)
                end = len(seg)+start
                mutant = self.prot[start:end]
                specificMutations[seg] = mutant
            else:
                mutationsPresent[seg]="SNP not in sequence range"
                
        self.mutationsPresent = mutationsPresent
        self.specificMutations = specificMutations
                

    def translateIndividualSeq(self, sequence, start):
        '''Takes a string as sequence instead of a list and translates only one 
        possible reading frame to amino acid sequence.  
        sequence: STRING, nucleic acid sequence
        start: INT, position to begin translation
        RETURN: STRING, translated amino acid sequence.
        '''
        codons = {'atg': 'M', 'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', 'ctt': 'L',
'ctc': 'L', 'cta': 'L', 'ctg': 'L','att': 'I', 'atc': 'I', 'ata': 'I', 'gtt': 'V',
'gtc': 'V', 'gta': 'V', 'gtg': 'V', 'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S',
'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', 'act': 'T', 'acc': 'T', 'aca': 'T',
'acg': 'T', 'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A', 'tat': 'Y', 'tac': 'Y', 
'taa': '*', 'tag': '*', 'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q', 'aat': 'N',
'aac': 'N', 'aaa': 'K', 'aag': 'K', 'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E',
'tgt': 'C', 'tgc': 'C', 'tga': '*', 'tgg': 'W', 'cgt': 'R', 'cgc': 'R', 'cga': 'R',
'cgg': 'R', 'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R', 'ggt': 'G', 'ggc': 'G',
'gga': 'G', 'ggg': 'G'} #change "STOP" to an asterisk
        wobble = start + 3   
        codonList = []
        codon = None
        translate = ''
        while wobble < len(sequence):
            codon = sequence[start:wobble]
            codonList.append(codon)
            try:
                translate+=codons[codon.lower()]
            except:
                translate+="N"
            start +=3
            wobble +=3
        return translate  
    
    def translator(self):
        '''Takes a sequence string and translates both forward and reverse possible amino
        acid sequences and output a dictionary of possible translations.
        INPUT: Sequence: STRING of sequenced nucleotides
        RETURN: DICTIONARY of possible translations.
        '''
        seq = Seq(self.sequence)
        revSeq = str(seq.reverse_complement())
        starts = [0, 1, 2]
        possibleTranslations = {}
        num = 1
        for start in starts:
            name = "+" + str(num)
            revName = "-" + str(num)
            translation = self.translateIndividualSeq(self.sequence, start)
            revTranslation = self.translateIndividualSeq(revSeq, start)
            possibleTranslations[name] = translation
            possibleTranslations[revName] = revTranslation
            num += 1
        return possibleTranslations
        
def mutationFinder(seq, targ):
    seq = resistanceMutationSearch(seq, targ)
    seq.determineTarget()
    
    seq.idCorrectProtein()
    if seq.prot != None:
        print seq.prot
        seq.detectMutations()
        print seq.mutationsPresent
        print seq.specificMutations
    else:
        print "Unable to translate sequence. Make sure sequence quality is good \n \
        and check to see if you entered the correct target."
        print "Target site selected: " + seq.targetSite

if __name__ == '__main__':
    mutationFinder(str(inFile.seq), target)