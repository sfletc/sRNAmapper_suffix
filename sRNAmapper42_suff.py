#     Copyright 2013 Stephen Fletcher Licensed under the
#     Educational Community License, Version 2.0 (the "License"); you may
#     not use this file except in compliance with the License. You may
#     obtain a copy of the License at
#
#      http://www.osedu.org/licenses/ECL-2.0
#
#     Unless required by applicable law or agreed to in writing,
#     software distributed under the License is distributed on an "AS IS"
#     BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
#     or implied. See the License for the specific language governing
#     permissions and limitations under the License.

import matplotlib
import numpy
from pylab import *
import matplotlib.pyplot as plt
import collections
from collections import Counter
import time
import argparse
import csv
from suffix_tree import *
      
class User_Input(object):
    """The command line class. Sets up the command line interface.
    """
    
    def __inti__(self):
        self._run = run
    
    def comline(self):
        """
        Allow script to be run from the command line or the GUI switched on
    
        Assumptions: all flags must be filled correctly or error
    
        comline() --> list(str, str, str, int, str, int, dtr ,int, int, str)
        """
    
                
        parser = argparse.ArgumentParser("python sRNAmapper42.py")
        
        
        parser.add_argument("-ana", "--analysis", type=str,
                help = "analysis type: dm, mntdm, mntcdp, cdp, cdps, cdpz, msfdm or msfcdp")
        
        parser.add_argument("-s1", "--seqfile1", type=str,
                help = "sequence file name 1")
                
        parser.add_argument("-s2", "--seqfile2", type=str,
                help = "sequence file name 2")

        parser.add_argument("-s3", "--seqfile3", type=str,
                help = "sequence file name 3")
        
        parser.add_argument("-s4", "--seqfile4", type=str,
                help = "sequence file name 4")
	                                           
        parser.add_argument("-r", "--reffile", type=str,
                help = "reference file name")
                
        parser.add_argument("-nt1", "--sRNAlength1", type=int,
                help = "sRNA nt 1 length")
                
        parser.add_argument("-nt2", "--sRNAlength2", type=int,
                help = "sRNA nt 2 length")
    
        parser.add_argument("-w", "--window_size", type=int,
                help = "window size for den_map_plot") 

        parser.add_argument("-a", "--absolute_den_map_plot", type=str,
                help = "y for absolute_den_map_plot") 

        parser.add_argument("-o", "--outfile", type=str,
                help = "outfilename")

        args = parser.parse_args()

        
        print "\n"

        print "Analysis type = "+str(args.analysis)
        print "\n"
        print "Reference File = "+ str(args.reffile)
        print "SeqFile 1 = "+ str(args.seqfile1)
        print "SeqFile 2 = "+ str(args.seqfile2)
        print "SeqFile 3 = "+ str(args.seqfile3)
        print "SeqFile 4 = "+ str(args.seqfile4)
        print "sRNA nt 1 length = "+ str(args.sRNAlength1)
        print "sRNA nt 2 length = "+ str(args.sRNAlength2)

        print "Output File = "+str(args.outfile)

        return [str(args.analysis), str(args.seqfile1), 
            'filler', str(args.seqfile2), 
            'filler', str(args.reffile),
            str(args.sRNAlength1), str(args.sRNAlength2), 
            str(args.outfile), 'filler', str(args.seqfile3), 
            'filler', str(args.seqfile4), 
            'filler','filler', 
            'filler','filler', 
            'filler','filler', 
            'filler','filler', 
            'filler','filler', 
            'filler','filler', 
            'filler','filler', 
            'filler','filler', 
            'filler', str(args.window_size), str(args.absolute_den_map_plot)]
                           

    def runner(self):
        """
        Run the selected analysis from the command line or switch to the GUI
        """
    
        run = self.comline()

    
        if run[0] == "dm":

            print "\nYou have selected the density map analysis\n"
           
            ana=Analysis()
            try:
                ana.dens_map(run[1], run[5], int(run[6]), run[8], str(run[31]), int(run[30]))
            except ValueError:
                ana.dens_map(run[1], run[5], int(run[6]), run[8], str(run[31]))
            
        elif run[0] == "mntdm":

            print "\nYou have selected the multi-sRNA length density map analysis\n"

            ana=Analysis()
            try:
                ana.multiplot_nt_dens_map(run[1], run[5], int(run[6]), int(run[7]),int(run[30]))
            except ValueError:
                ana.multiplot_nt_dens_map(run[1], run[5], int(run[6]), int(run[7])) 
        elif run[0] == "msfdm":
            
            print "\nYou have selected the multi-sequence file density map analysis\n"
            
            ana=Analysis()
            try:
                ana.multiplot_seq_dens_map(run[1], run[2], run[3], run[4], run[5], 
                    int(run[6]),int(run[30]))
            except ValueError:
                ana.multiplot_seq_dens_map(run[1], run[2], run[3], run[4], run[5], 
                    int(run[6]))
        elif run[0] == "cdp":
            
            print "\nYou have selected the comparitive density plot analysis\n"
            
        
            ana = Analysis()
            ana.comp_reads(run[1], run[3], run[5], int(run[6]), 
                run[8])
        
        elif run[0] == "cdps":
            
            print "\nYou have selected the comparitive density plot analysis with sRNA reference\n"
            
        
            ana = Analysis()
            ana.comp_reads_sRNA(run[1], run[2], run[3], run[4], run[5], 
                run[8])

        elif run[0] == "cdpz":
            
            print "\nYou have selected the comparitive density plot inc. zero analysis\n"
            
        
            ana = Analysis()
            ana.comp_reads_zero(run[1], run[2], run[3], run[4], run[5], int(run[6]), 
                run[8])


        elif run[0] =="msfcdp":
            
            print "\nYou have selected the multi-sequence file comparitive density plot analysis\n"
            
            
            ana=Analysis()
            ana.comp_seqfile_multi(run[1], run[3], run[10], run[12], run[5], int(run[6]))
        
        elif run[0] =="mntcdp":
            
            print "\nYou have selected the multi-nucleotide comparitive density plot analysis\n"
            
            ana=Analysis()
            ana.comp_reads_2xnt(run[1], run[3], run[5], 
                int(run[6]), int(run[7]), run[8])   

class Analysis(object):
    """For running the analyses
    """
    
    def __init__(self):
        self._ref_dict1 = Ref_Dict()
        self._ref_dict2 = Ref_Dict()
        self._seq_dict1 = Seq_Dict()
        self._seq_dict2 = Seq_Dict()
        self._seq_dict3 = Seq_Dict()
        self._seq_dict4 = Seq_Dict()   
        self._seq_dict1a = Seq_Dict()
        self._seq_dict1b = Seq_Dict()
        self._seq_dict2a = Seq_Dict()
        self._seq_dict2b = Seq_Dict()
        self._samp1=Align_Seq()
        self._samp2=Align_Seq()
        self._samp3=Align_Seq()
        self._samp4=Align_Seq()
        self._samp1a=Align_Seq()
        self._samp1b=Align_Seq()
        self._samp2a=Align_Seq()
        self._samp2b=Align_Seq()
        
    
    def dens_map(self, seq1, ref1, nt1, out_file, absol='n', w=50):
        """
        Analysis that plots a density map as counts at each alignment position
        on a single reference sequence, positive for the forward strand and
        negative for the reverse strand
        
        dens_map(str, str, int, str) --> None
        """
        a='n'
        if absol is None: 
        	a='n'
    	elif absol =='y':
    		a='y'
        if seq1 =='' or ref1=='' or nt1==0 or out_file=='':
            raise Incomplete_Args
            
        seq_dict= self._seq_dict1.get_seq(seq1, nt1)
        ref_dict= self._ref_dict1.get_ref(ref1)
        start = time.clock()
        samp1 = self._samp1.align_seq(seq_dict[0], ref_dict, nt1).hit_counter()
        out_samp1 = Out_Put(samp1)
        out_samp1.write2csv(out_file).den_map_plot(seq1, ref1, a, w)  
        

        print ".........we're done\n"
    
        print "Time taken = " + str((time.clock() - start)) +" seconds"


    
    def multiplot_nt_dens_map(self, seq1, ref1, nt1, nt2, w=50):
        """
        Analysis that plots a density map as counts at each alignment position
        on a single reference sequence, positive for the forward strand and
        negative for the reverse strand and in different colours for each sRNA
        length
        
        dens_map(str, str, int, int, str) --> None
        """
        
        if seq1 =='' or ref1=='' or nt1==0 or nt2 ==0:
            raise Incomplete_Args
                    
        seq_dict1= self._seq_dict1.get_seq(seq1, nt1)
        seq_dict2= self._seq_dict1.get_seq(seq1, nt2)
        ref_dict1= self._ref_dict1.get_ref(ref1)
        
        start = time.clock()
        
        samp1 = self._samp1.align_seq(seq_dict1[0], ref_dict1, nt1).hit_counter() 
        samp2 = self._samp2.align_seq(seq_dict2[0], ref_dict1, nt2).hit_counter()   
        
        out_samp1=Out_Put(samp1)
        out_samp2=Out_Put(samp2)
        
        
        out_samp1.den_multi_plot(out_samp2, (str(nt1)+" nucleotides"), 
            (str(nt2)+" nucleotides"), ref1, "sRNA Density", w)
        print ".........we're done\n"
        
        
        print "Time taken = " + str((time.clock() - start)) +" seconds"


    def multiplot_seq_dens_map(self, seq1, count1, seq2, count2, ref1, 
        nt1, w=50):
        """
        Analysis that plots a density map as counts at each alignment position
        on a single reference sequence, positive for the forward strand and
        negative for the reverse strand and in different colours for each 
        sequence file
        
        dens_map(str, int, str, int, str, int, str) --> None
        """
        
        if seq1 =='' or count1==0 or seq2 == '' or ref1=='' or nt1==0:
            raise Incomplete_Args      
                          
        if seq1 =='' or count1 == 0 or seq2 =='' or count2 == 0 \
            or ref1=='' or nt1==0:
            raise Incomplete_Args   
            
        
        start = time.clock()
        
        seq_dict1= self._seq_dict1.get_seq(seq1, nt1)
        seq_dict2= self._seq_dict2.get_seq(seq2, nt1)
        ref_dict1= self._ref_dict1.get_ref(ref1)
        
        samp1 = self._samp1.align_seq(seq_dict1[0], ref_dict1, nt1).\
            hit_counter().list_norm(seq_dict1[1])
        samp2 = self._samp2.align_seq(seq_dict2[0], ref_dict1, nt1).\
            hit_counter().list_norm(seq_dict2[1])  
        

        out_samp1=Out_Put(samp1)
        out_samp2=Out_Put(samp2)
        
        out_samp1.den_multi_plot(out_samp2, seq1, seq2, ref1, 
            "sRNA counts per million total counts")
    
        print ".........we're done\n"
        
    
        print "Time taken = " + str((time.clock() - start)) +" seconds"        
        

    def comp_reads(self, seq1, seq2, ref1, nt1, 
        outfile):
        """
        Analysis that plots the counts of alignments from two sequence files 
        vs. sequences in one reference file as x,y co-ordinates.
        
        comp_reads(str, int, str, int, str, int, str) --> None
        """
                
        if seq1 =='' or seq2 =='' \
            or ref1=='' or nt1==0 or outfile=='':
            raise Incomplete_Args    
        
        start = time.clock()
        
        ref_dict_1 = self._ref_dict1.get_ref(ref1)
        seq_dict_1 = self._seq_dict1.get_seq(seq1, nt1)
        seq_dict_2 = self._seq_dict2.get_seq(seq2, nt1)
       
        samp1 = self._samp1.align_seq(seq_dict_1[0], 
            ref_dict_1, nt1).hit_counter().hit_per_seq().normit(seq_dict_1[1])   
        
        print ".........spinning up the hamsters"
        samp2 = self._samp2.align_seq(seq_dict_2[0], ref_dict_1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict_1[1])
   
        print ".........launching the llamas"
    
        comp_both = samp1.comp_hit_per_seq(samp2)
        comp_both.write_coords_csv(outfile)
        comp_both.cdp_plotter(seq1,seq2)
         
        print "Time taken = " + str((time.clock() - start)) +" seconds"

    def comp_reads_sRNA(self, seq1, seq1_total, seq2, seq2_total, ref1, 
        outfile):
        """
        Analysis that plots the counts of alignments from two sequence files 
        vs. sequences in one reference file as x,y co-ordinates.
        
        comp_reads(str, int, str, int, str, int, str) --> None
        """
                
        if seq1 =='' or seq1_total == 0 or seq2 =='' or seq2_total == 0 \
            or ref1==''  or outfile=='':
            raise Incomplete_Args    
        
        start = time.clock()
        
        ref_dict_1 = self._ref_dict1.get_ref_sRNA(ref1)
        seq_dict_1 = self._seq_dict1.get_seq_sRNA(seq1)
        seq_dict_2 = self._seq_dict2.get_seq_sRNA(seq2)
        
        samp1 = self._samp1.align_seq_sRNA(seq_dict_1, 
            ref_dict_1).hit_per_seq().normit(seq1_total)   
        
        print ".........spinning up the hamsters"
        samp2 = self._samp2.align_seq_sRNA(seq_dict_2, 
            ref_dict_1).hit_per_seq().normit(seq2_total)
   
        print ".........launching the llamas"
    
        comp_both = samp1.comp_hit_per_seq(samp2)
        comp_both.write_coords_csv(outfile)
        comp_both.cdp_plotter(seq1,seq2)
        
   
        print ".........we're done\n"
        print 
    
        print "Time taken = " + str((time.clock() - start)) +" seconds"


    # def comp_reads_zero(self, seq1, seq1_total, seq2, seq2_total, ref1, nt1, 
    #     outfile):
    #     """
    #     Analysis that plots the counts of alignments from two sequence files 
    #     vs. sequences in one reference file as x,y co-ordinates.
        
    #     comp_reads(str, int, str, int, str, int, str) --> None
    #     """
                
    #     if seq1 =='' or seq1_total == 0 or seq2 =='' or seq2_total == 0 \
    #         or ref1=='' or nt1==0 or outfile=='':
    #         raise Incomplete_Args    
        
    #     start = time.clock()
        
    #     ref_dict_1 = self._ref_dict1.get_ref(ref1)
    #     seq_dict_1 = self._seq_dict1.get_seq(seq1, nt1)
    #     seq_dict_2 = self._seq_dict2.get_seq(seq2, nt1)
       
    #     samp1 = self._samp1.align_seq(seq_dict_1, 
    #         ref_dict_1, nt1).hit_counter().hit_per_seq().normit(seq1_total)   
        
    #     print ".........spinning up the hamsters"
    #     samp2 = self._samp2.align_seq(seq_dict_2, ref_dict_1, 
    #         nt1).hit_counter().hit_per_seq().normit(seq2_total)
   
    #     print ".........launching the llamas"
    
    #     comp_both = samp1.comp_hit_per_seq_zero(samp2)
    #     comp_both.write_coords_csv(outfile)
    #     comp_both.cdp_plotter(seq1,seq2)
        
   
    #     print ".........we're done\n"
    #     print 
    
    #     print "Time taken = " + str((time.clock() - start)) +" seconds"


    def comp_seqfile_multi(self, seq1, seq2, seq3, seq4, ref1, nt1):
        
        start = time.clock()

        seq_dict1 = self._seq_dict1.get_seq(seq1, nt1)
        seq_dict2 = self._seq_dict2.get_seq(seq2, nt1)
        seq_dict3 = self._seq_dict3.get_seq(seq3, nt1)
        seq_dict4 = self._seq_dict4.get_seq(seq4, nt1)
        ref_dict1 = self._ref_dict1.get_ref(ref1)

        samp_seq1 = self._samp1.align_seq(seq_dict1[0], ref_dict1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict1[1])
        self._samp1=Align_Seq()

        samp_seq2 = self._samp2.align_seq(seq_dict2[0], ref_dict1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict2[1])
        self._samp2=Align_Seq()

        samp_seq3 = self._samp3.align_seq(seq_dict3[0], ref_dict1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict3[1])
        self._samp1=Align_Seq()

        samp_seq4 = self._samp4.align_seq(seq_dict4[0], ref_dict1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict4[1])
        self._samp1=Align_Seq()

        comp_a = samp_seq1.comp_hit_per_seq(samp_seq2)
        comp_b = samp_seq3.comp_hit_per_seq(samp_seq4)

        comp_a.cdp_2seq_plotter(comp_b, seq1, seq2, seq3, seq4)

        print "\nTime taken = " + str((time.clock() - start)) +" seconds"
        
    
    def comp_reads_2xnt(self, seq1, seq2, ref1, nt1, 
        nt2, out_file):
        """
        Analysis that plots the counts of alignments from two sequence files 
        vs. sequences in one reference file as x,y co-ordinates for two sRNA
        lengths in different colours.
        
        comp_reads_multi(str, int, str, int, str, int, int, str) --> None
        """
        
        if seq1 =='' or ref1=='' or nt1==0 or out_file=='':
            raise Incomplete_Args
        
        start = time.clock()
        
        seq_dict1_nt1= self._seq_dict1a.get_seq(seq1, nt1)
        seq_dict2_nt1= self._seq_dict1b.get_seq(seq2, nt1)
        seq_dict1_nt2= self._seq_dict2a.get_seq(seq1, nt2)
        seq_dict2_nt2= self._seq_dict2b.get_seq(seq2, nt2)
        ref_dict_1= self._ref_dict1.get_ref(ref1)
        
        samp_seq1_nt1 = self._samp1a.align_seq(seq_dict1_nt1[0], ref_dict_1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict1_nt1[1])
        self._samp1a=Align_Seq()   
        
        samp_seq2_nt1 = self._samp1b.align_seq(seq_dict2_nt1[0], ref_dict_1, 
            nt1).hit_counter().hit_per_seq().normit(seq_dict2_nt1[1])
        self._samp1b=Align_Seq()    
        
        samp_seq1_nt2 = self._samp1a.align_seq(seq_dict1_nt2[0], ref_dict_1, 
            nt2).hit_counter().hit_per_seq().normit(seq_dict1_nt2[1])
        samp_seq2_nt2 = self._samp1b.align_seq(seq_dict2_nt2[0], ref_dict_1, 
            nt2).hit_counter().hit_per_seq().normit(seq_dict2_nt2[1])
    
        comp_both_2xnt_list = samp_seq1_nt1.mutli_comp_hit_per_seq_2xnt\
            (samp_seq2_nt1, samp_seq1_nt2, samp_seq2_nt2)
       
        (Out_Put(comp_both_2xnt_list[0])).cdp_2nt_plotter\
            (Out_Put(comp_both_2xnt_list[1]), seq1, seq2, nt1, nt2) 
        
        print "\nTime taken = " + str((time.clock() - start)) +" seconds"     
                   

                        

class Dicter(dict):
    """
    This class is set up to as the main class for the reference and sequence
    dictionaries
    """
    
    def __init__(self):
        self._sequence = ''
        self._header = ''
        

class Ref_Dict(Dicter):
    """The reference dictionary class stores reference information as a header
    and sequence.  Each individual reference has a unique header which, in fasta
    format, is always in indicated by the symbol >
    """
    
    def __init__(self):
        Dicter.__init__(self)
        self.ref_dict = {}
        
    
    def get_ref(self, ref_filename):
        """Load file in fasta format.  Make a dictionary with header as key and
        sequence as value.  Initially, f if added to the beginning
        of the header (key) to indicate the sequence is in the 5' to 3' 
        direction. The complement is then calculated (opposite strand 5' to 3')
        and an r added to its header (key).
    
        Assumptions: ref_filename is a text file in fasta format.  Sequences
        are DNA.
    
        get_ref(str) --> None
        """
        
        self.filename = Dict_File(ref_filename)
        self.loaded_ref = self.filename.load_dict()
                
        full_len_seq = ''
        key=''
        for line in self.loaded_ref:
            
            clean_line=line.strip()
            
            if line[0] == '>' and full_len_seq == '':
                key = clean_line
            elif line[0] == '>' and full_len_seq != '':
                self.ref_dict.update({'f'+ key:DNA(full_len_seq)})
                self.ref_dict.update({'r'+ key:(DNA(full_len_seq)).
                    complement()})
                key=clean_line
                full_len_seq = ''
            
            elif line[0] == '' and full_len_seq != '':
                self.ref_dict.update({'f'+key:DNA(full_len_seq)})
                self.ref_dict.update({'r'+key:(DNA(full_len_seq)).complement()})
                key=clean_line
                full_len_seq = ''   
            
            elif line[0] =='':
                pass
            
            else: 
                full_len_seq += clean_line.upper()
        
        self.ref_dict.update({'f'+key:DNA(full_len_seq)})
        self.ref_dict.update({'r'+key:(DNA(full_len_seq)).complement()}) 
        
        
        return self.ref_dict
        
    def suffix_ref(self, ref_filename):
        self.filename = Dict_File(ref_filename)
        self.loaded_ref = self.filename.load_dict()
                
        full_len_seq = ''
        key=''
        for line in self.loaded_ref:
            
            clean_line=line.strip()
            
            if line[0] == '>' and full_len_seq == '':
                key = clean_line
            elif line[0] == '>' and full_len_seq != '':
                self.ref_dict.update({'f'+ key:SuffixTree(DNA(full_len_seq))})
                self.ref_dict.update({'r'+ key:SuffixTree(DNA(full_len_seq).complement())})
                key=clean_line
                full_len_seq = ''
            
            elif line[0] == '' and full_len_seq != '':
                self.ref_dict.update({'f'+ key:SuffixTree(DNA(full_len_seq))})
                self.ref_dict.update({'r'+ key:SuffixTree(DNA(full_len_seq).complement())})
                key=clean_line
                full_len_seq = ''   
            
            elif line[0] =='':
                pass
            
            else: 
                full_len_seq += clean_line.upper()
        
        self.ref_dict.update({'f'+ key:SuffixTree(DNA(full_len_seq))})
        self.ref_dict.update({'r'+ key:SuffixTree(DNA(full_len_seq).complement())})
        return self.ref_dict

    def get_ref_sRNA(self, ref_filename):
        """Load file in fasta format.  Converts U to T (ie RNA to DNA).
        Make a dictionary with header as key and
        sequence as value.  Initially, f if added to the beginning
        of the header (key) to indicate the sequence is in the 5' to 3' 
        direction. No complement calculated.
    
        Assumptions: ref_filename is a text file in fasta format.  Sequences
        are DNA.
    
        get_ref(str) --> None
        """

        self.filename = Dict_File(ref_filename)
        self.loaded_ref = self.filename.load_dict()
                
        full_len_seq = ''
        key=''
        for line in self.loaded_ref:
            
            clean_line=line.strip()
            
            if line[0] == '>' and full_len_seq == '':
                key = clean_line
            elif line[0] == '>' and full_len_seq != '' and DNA(full_len_seq) not in self.ref_dict.values():
                self.ref_dict.update({'f'+ key:DNA(full_len_seq)})
                # self.ref_dict.update({'r'+ key:(DNA(full_len_seq)).
                #     complement()})
                key=clean_line
                full_len_seq = ''
            
            elif line[0] == '' and full_len_seq != '' and DNA(full_len_seq) not in self.ref_dict.values():
                self.ref_dict.update({'f'+key:DNA(full_len_seq)})
                #self.ref_dict.update({'r'+key:(DNA(full_len_seq)).complement()})
                key=clean_line
                full_len_seq = ''   
            
            elif line[0] =='':
                pass
            
            else: 
                
                full_len_seq=clean_line.upper().replace("U","T")
            
        
        self.ref_dict.update({'f'+key:DNA(full_len_seq)})
        #self.ref_dict.update({'r'+key:(DNA(full_len_seq)).complement()}) 

        return self.ref_dict

    def __repr__(self):
        return "Ref_Dict_({0})".format(self.ref_dict)
        
    def __str__(self):
        return "Ref_Dict: {0}".format(self.ref_dict)


class Seq_Dict(Dicter):
    """
    The Sequence Dictionary class stores the sequence and count of the sRNA. 
    The sequence is unique is forms the key and the count the value.  The sRNA 
    is always in DNA format (this is how sRNAs are sequenced).  Only sRNAs
    of the selected length are added to the dictionary.
    """
    
    def __init__(self):
        Dicter.__init__(self)
        self.seq_dict = {}
    

    def get_seq(self, seq_filename, nt):
        """
        Loads the sequence file. Only includes sequences of length nt.
        
        Asumption: text file is in BGI format: sRNA length, sRNA count, sRNA seq

        get_seq(str, int) --> None
        """
        
        self.filename = Dict_File(seq_filename)
        self.loaded_seq = self.filename.load_dict()
        count = 0
        for i in self.loaded_seq:
        
            a = i.strip() 

            b = a.split('\t')
            count += int(b[1])
            if int(b[0]) == nt:
                seq = DNA(b[2])
                self.seq_dict.update({seq:(b[1])})
        
        return self.seq_dict, count

    def get_seq_sRNA(self, seq_filename):
        """
        

        Loads the sequence file. Only includes sequences of length nt.
        
        Asumption: text file is in BGI format: sRNA length, sRNA count, sRNA seq

        get_seq(str, int) --> None
        """
        
        self.filename = Dict_File(seq_filename)
        self.loaded_seq = self.filename.load_dict()

  
        for i in self.loaded_seq:
        
            a = i.strip() 

            b = a.split('\t')

            
            seq = DNA(b[2])
            self.seq_dict.update({seq:(b[1])})
        
        return self.seq_dict
                                                                                          
class Dict_File(Dicter):
    """
    Sets up a Dict_File in the analysis pipeline
    """
    
    def __init__(self, file_name):
         self.file_name = file_name
         
    def load_dict(self):
         self.loaded_ref = open(self.file_name, 'rU')
         
         return self.loaded_ref
         
     
    def __repr__(self):
         return "Dict_File({0})".format(self.file_name)
        
    def __str__(self):
        return "Dict_File: {0}".format(self.file_name)
                        

class DNA(str):
    """
    A DNA sequence is essentially a series of A, G, C and Ts.  Thus, a subclass
    of string seems appropriate.
    """
    
    def __init__(self, sequence):
        self._sequence = sequence
        self._header = ''
        self._nuc = ('A', 'G', 'C', 'T')
        
          
    def add_header(self, header):
        self._header = header

    def get_header(self):
        return self._header   
    
    def check_all_DNA(self):
        for i in self._sequence:
            if i not in self._nuc:
                return False
            else:
                return True
     
    def complement(self):
        """Provides the complement in the 5' - 3' direction

        Assumption: reference consists of A, G, C, T only

        complement(str) --> str
        """
        d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return DNA(''.join(d[c] if c in d else c for c in reversed(self)))
        
    
    def __len__(self):
        return len(self._sequence)
       
            
    def __repr__(self):
        return "DNA({0})".format(self._sequence)
        
    def __str__(self):
        return "DNA sequence: {0}".format(self._sequence)
    
    def __getslice__(self, i, j):
        #return DNA(self[max(0, i):max(0, j):]) 
        return DNA(self[i:j:]) 

class Align_Seq(object):
    """
    Align_Seq is the central class for alignment.  An object of the Align_Seq
    class can be acted on by the align_seq method.      
    """

    def __init__(self):
        self._results = []
        self._ordered = []
             
    def align_seq(self, seq_dict1, ref_dict1, nt):
        """
        The align seq method checks if a sequence of a set length from the
        seq_dict is the same as a sequence in a sequence in the ref_dict by
        scanning the ref_dict sequence in a sliding window of that length.
        
        Returned is a sorted list of header, sRNA, count, nt length 
        and direction
        
        Assumptions: seq and ref dicts are in the correct format.
        
        align_seq(dict, dict, int) --> Post_Process(list
                                                (str, str, int, int, str))
        """

        for header, seq in ref_dict1.items():
      
            count_start = 0
        
            while count_start < (len(seq)-(nt-1)):
                   
                query_seq = seq[count_start:(count_start+nt)]
                if query_seq in seq_dict1 and header[0] =='f':
                    self._results.append([header, query_seq, 
                        seq_dict1[query_seq], count_start, 'f'])

                elif query_seq in seq_dict1 and header[0] =='r':
                    self._results.append([header, query_seq, 
                        seq_dict1[query_seq], len(seq)-count_start, 'r']) #check if start point correct

                else:
                    pass
                        
                count_start +=1
                                
        self._ordered = sorted(self._results, key=lambda _results:_results[1])
        return Post_Process(self._ordered)

    def align_seq_sRNA(self, seq_dict1, ref_dict1):
        """need to remove nt from calling functions as not needed.  Also f and count_start
        """

        for header, seq in ref_dict1.items():

            if seq in seq_dict1:
                self._results.append([header, seq, int(seq_dict1[seq]), 0, 'f'])

            else:
                pass
        self._ordered = sorted(self._results, key=lambda _results:_results[1])

        return Post_Process(self._ordered)

class Post_Process(Align_Seq):
    """Objects in the Post_Process class are processed in a variety of ways
    """
    
    def __init__(self, ordered_list):
        Align_Seq.__init__(self)
        self._ordered_list = ordered_list
    
    def __iter__(self):
        return iter(self._ordered_list)
    
    def hit_counter(self):
        """
        If an sRNA aligns to multiple reference sequences, divides the count per
        alignment by the number of sequences aligned to.  ie. an sRNA with a 
        count of 50 that aligns to 2 ref seqs has 25 assigned to each.
        """
        count=0
        _sRNA_mapped_count = 0

        _ordered_list = self._ordered_list
        hit_counter = collections.Counter()
    
        while count < len(_ordered_list):
            _ordered_list[count][0]=_ordered_list[count][0][1:]
            hit_counter.update(_ordered_list[count][1:2]) #query seq and count
            
            count+=1
            
        hit_counter = dict(hit_counter) #makes it a dictionary 
        
    
        for i in _ordered_list:
            
            i[2] = float(i[2])/float(hit_counter[i[1]]) 
                #divides the count for an alignment by the alignments  
            _sRNA_mapped_count += i[2]
        self._ordered = sorted(_ordered_list, 
            key=lambda _ordered_list:_ordered_list[0])
        print str(int(_sRNA_mapped_count)) +" reads have mapped"
        return Post_Process(self._ordered)


    def list_norm(self, total_cnt):
        """
        Normalises the alignment count to hits per million by using the 
        total number of reads per sequence file.  Only needed when analysis
        involves two sequence files, as ensures the results can be directly 
        compared.  
        
        list_norm(list(str, str, int, int, str), int) --> None
        """
        
        total_cnt = float(total_cnt)
        per_million = total_cnt/1000000
        
        for i in self._ordered_list:
            i[2] = i[2]/per_million
        
        return self._ordered_list
            
    
    def hit_per_seq(self):
        """
        Returns total number of all alignments per reference seq.
        """
        
        
        total = 0
        self._dicter = {self._ordered_list[0][0]:self._ordered_list[0][2]}
        for a in self._ordered_list:
            if a[0] in self._dicter:
                total += a[2]
                self._dicter[a[0]] = total
                
            else:
                total = a[2]
                self._dicter[a[0]] = total
                
        return Normit(self._dicter)


class Normit(Post_Process):
    """Normalisation class for dictionaries
    """
    
    def __init__(self, dicter):
        self._dicter = dicter    
    
    def normit(self, total_cnt):
        """
        Normalises the alignment count to hits per million by using the 
        total number of reads per sequence file.  Only needed when analysis
        involves two sequence files, as ensures the results can be directly 
        compared.  
        
        normit(dict(str:int), int) --> None
        """
        
        total_cnt = float(total_cnt)
        per_million = total_cnt/1000000
        
        for i in self._dicter:
            
            self._dicter[i] = (self._dicter[i]/per_million)
        return Comp_Hit_Per_Seq(self._dicter)


class Comp_Hit_Per_Seq(Post_Process):
    """
    Objects of this class are used for comparative  analyses
    """
    
    def __init__(self, dicter):
        self._dicter = dicter
    def __iter__(self):
        return iter(self._dicter)
    def __getitem__(self, k):
        return self._dicter[k]
    
    def comp_hit_per_seq(self, dicter2):
        """
        Takes two alignment dictionaries (ref seq:alignments) and returns counts
        as x,y co-ordinates in a tuple only when a reference sequences has 
        alignments in both dictionaries
        
        comp_hit_per_seq(dict(str:int), dict(str:int)) --> (int,int)
        """
        
        co_ord = []
        hits1_hits2 = ()
        for key in self._dicter:
            if key in dicter2:
                hits1_hits2 = (self._dicter[key], dicter2[key])
                co_ord.append(hits1_hits2)
        
        return Out_Put(co_ord)

    def comp_hit_per_seq_zero(self, dicter2):
        """
        Takes two alignment dictionaries (ref seq:alignments) and returns counts
        as x,y co-ordinates in a tuple even if alignment only in one dictionary
        
        comp_hit_per_seq(dict(str:int), dict(str:int)) --> (int,int)
        """
        
        co_ord = []
        hits1_hits2 = ()
        for key in self._dicter:
            if key in dicter2:
                hits1_hits2 = (self._dicter[key], dicter2[key])
                co_ord.append(hits1_hits2)
            else:
                hits1_hits2 = (self._dicter[key], 0)
                co_ord.append(hits1_hits2)
        hits1_hits2 = ()
        for key in dicter2:
            if key not in self._dicter:
                hits1_hits2 = (0, dicter2[key])
                co_ord.append(hits1_hits2)
        
        return Out_Put(co_ord)
    
    def multi_comp_hit_per_seq(self, samp_seq2_nt1, samp_seq1_nt2, 
        samp_seq2_nt2):
        """
        Takes four alignment dictionaries (ref seq:alignments) and returns 
        counts as x,y co-ordinates in a tuple only when a reference sequences 
        has alignments in all dictionaries and the count is not 0 in any
        dictionary
        
        comp_hit_per_seq(dict(str:int), dict(str:int), 
                dict(str:int), dict(str:int)) --> (int,int)
        """
        
        co_ord=[]
        hits1_hits2=()
        for key in self._dicter:

            if key in samp_seq2_nt1 and key in samp_seq1_nt2 and key in \
                samp_seq2_nt2 and self._dicter[key] != 0 and \
                samp_seq2_nt1[key]!=0 and samp_seq1_nt2[key] !=0 \
                and samp_seq2_nt2[key] !=0:
                hits1_hits2 = ((self._dicter[key] / samp_seq2_nt1[key]),
                (samp_seq1_nt2[key]/samp_seq2_nt2[key]))
                co_ord.append(hits1_hits2)
                
        return Out_Put(co_ord)       
    

    def mutli_comp_hit_per_seq_2xnt(self, samp_seq2_nt1, samp_seq1_nt2, 
        samp_seq2_nt2):
        """
        For multiple comparative analysis (ie. comparative density plot with two
        sRNA lengths) but want separate plotting. Takes four dictionaries and 
        returns two lots of x, y coordinates with a reference sequence has 
        alignments in all four dictionaries
        
        multi_comp_hit_per_seq_2xnt(dict(str:int), dict(str:int), 
                dict(str:int), dict(str:int)) --> [(int,int), (int,int)]
        """
        
        co_ord_nt1=[]
        co_ord_nt2=[]
        hits_all_nt1=()
        hits_all_nt2=()
        out_put_list=[]

        
        for key in self._dicter:
            
            
            if key in samp_seq2_nt1 and key in samp_seq1_nt2 and key in \
                samp_seq2_nt2:
                
                hits_all_nt1 = (self._dicter[key], samp_seq2_nt1[key])
                
                hits_all_nt2 = (samp_seq1_nt2[key], samp_seq2_nt2[key])
               
                co_ord_nt1.append(hits_all_nt1)
                
                co_ord_nt2.append(hits_all_nt2)
                
        out_put_list.append(co_ord_nt1)
        out_put_list.append(co_ord_nt2)
        
        return out_put_list
        

class Out_Put(object):
    """
    Objects in the Out_Put class are used for writing to .csv files and plotting
    """
    
    def __init__(self, results_list):
        self._results_list = results_list
    def __iter__(self):
        return iter(self._results_list)  
    
    
    def write2csv(self, outfile):
        """
        Write list to .csv for excel in delineated  format
    
        write2csv(list(str, str, int, int, str), str) 
            --> list(str, str, int, int, str)
        
        """
                
        mycsv = csv.writer(open(outfile, 'wb'))
        headers = ('Ref_seq_ID', 'sRNA_seq', 'count', "pos_from_5'_fwd",
                                                            "Fwd_or_Rvs")
        mycsv.writerow(headers)
        
        for row in self._results_list:
            if row[4] == 'r': row[2] = -row[2]
            mycsv.writerow(row)
            if row[4] == 'r': row[2] = -row[2]
        return Out_Put(self._results_list)

    
    def write_coords_csv(self, outfile):
        """
        Writes list of tuples to excel in x,y columns
        
        write_coords_csv(list(int, int), str) --> list(int,int)
        """
            
        mycsv = csv.writer(open(outfile, 'wb'))
        headers = ('Sample 1 hits', 'Sample 2 hits')
        mycsv.writerow(headers)
        for i in self._results_list:
            row = [i[0], i[1]]
            mycsv.writerow(row)  
        return Out_Put(self._results_list)
    
    def cdp_plotter(self, seq1, seq2):
        """
        Plots a comparative  density plot (ie. 1 pixel dot at x,y)
        
        cdp_plotter(list(x,y), str, str) --> None
        Interactive mode so displays on Mac OS X properly - remove for Win
        """
        self._results_list=sorted(self._results_list)
     
        _max=max(self._results_list[-1][0],self._results_list[-1][1])
        
        plt.scatter(*zip(*self._results_list), s=3)
        
        arrow(0.1,0.1,5000000,5000000, color = 'r')
        xlabel(seq1)
        ylabel(seq2)
        xscale('log')
        yscale('log')
        xlim(1,_max+100)
        ylim(1,_max+100)
        plt.title('Comparative Density Plot')
        plt.ion() #Remove if running on windows
        plt.show()
        
     
    def cdp_2nt_plotter(self, ntx2_result, seq1, seq2, nt1, nt2):
        """
        Plots a multi-comparative  density plot (ie. 1 pixel dot at x,y) in 
        two colours representing each sRNA length
        
        cdp_2nt_plotter(list(x,y), str, str, int, int) --> None
        Interactive mode so displays on Mac OS X properly - remove for Win
        """
        self._results_list=sorted(self._results_list)
        ntx2_result=sorted(ntx2_result)

        _max=max(self._results_list[-1][0],self._results_list[-1][1], ntx2_result[-1][0], ntx2_result[-1][1])

        leg_1 = str(nt1)+' nucleotides'
        leg_2 = str(nt2)+' nucleotides'
        plt.scatter(*zip(*self._results_list), s=1, color = 'r')
        plt.scatter(*zip(*ntx2_result), s=1, color = 'g')
        arrow(1,1,10000,10000, color = 'r')
        xlabel(seq1)
        ylabel(seq2)
        xscale('log')
        yscale('log')
        xlim(1,_max+100)
        ylim(1,_max+100)
        legend([leg_1,leg_2],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
        plt.title('Multi-nucleotide comparative density plot')
        plt.ion() #Remove if running on windows
        plt.show()
    
    def cdp_2seq_plotter(self, result_2, seq1, seq2, seq3, seq4):
        """
        Plots a multi-comparative  density plot (ie. 1 pixel dot at x,y) in 
        two colours representing each sRNA length
        
        cdp_2nt_plotter(list(x,y), str, str, int, int) --> None
        Interactive mode so displays on Mac OS X properly - remove for Win
        """
        self._results_list=sorted(self._results_list)
        result_2=sorted(result_2)
        _max=max(self._results_list[-1][0],self._results_list[-1][1], result_2[-1][0], result_2[-1][1])
        leg_1 = str(seq1)+ ' ' +str(seq2)
        leg_2 = str(seq3)+ ' ' +str(seq4)
        plt.scatter(*zip(*self._results_list), s=1, color = 'r')
        plt.scatter(*zip(*result_2), s=1, color = 'g')
        arrow(1,1,10000,10000, color = 'r')
        xlabel(seq1 + ' + ' + seq3)
        ylabel(seq2 + ' + ' + seq4)
        xscale('log')
        yscale('log')
        xlim(1,_max+100)
        ylim(1,_max+100)
        legend([leg_1,leg_2],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
        plt.title('Multi-sequence file comparative density plot')
        plt.ion() #Remove if running on windows
        plt.show()   

    def org_coord(self):
        """
        For plotting - so can use smoothed window, returns x co-ords, y co-ords(fwd strand)
        and y co-ords (rev strand) as lists
        """
        plot_ord_list_f = []
        plot_ord_list_r = []
        for result in self._results_list:
            if result[4] == 'r':
                plot_ord_list_r.append((result[3],(0-result[2])))
            else:
                plot_ord_list_f.append((result[3],(result[2])))        
        plot_ord_list_f = sorted(plot_ord_list_f)
        plot_ord_list_r = sorted(plot_ord_list_r)

        pos_f=0
        pos_r=0
        i=0
        y_c_f=[]
        y_c_r=[]
        x_c=range(min(plot_ord_list_f[0][0], plot_ord_list_r[0][0]),max(1+plot_ord_list_f[-1][0],\
            1+plot_ord_list_r[-1][0]))
        
        while i < len(x_c):
                try:
                    if x_c[i]==plot_ord_list_f[i-pos_f][0]:
                            y_c_f.append(plot_ord_list_f[i-pos_f][1])
                    else:
                            y_c_f.append(0)
                            pos_f +=1
                except IndexError:
                    y_c_f.append(0)
                try:

                    if x_c[i]==plot_ord_list_r[i-pos_r][0]:
                            y_c_r.append(plot_ord_list_r[i-pos_r][1])
                            
                    else:
                            y_c_r.append(0)
                            pos_r +=1

                except IndexError:
                    y_c_r.append(0) #should be ycr probably
                i+=1

        return x_c, y_c_f,y_c_r

    def den_map_plot(self, seq1, ref1, abs='n',w=50):
        """
        Plots a density map plot
        
        den_map_plot(list(str, str, int, int, str), str, str) --> None
        """
        #abs='n'

        a=self.org_coord()
        x_c= a[0]
        y_c_f=a[1]
        y_c_r=a[2]    
        if abs=='n':

            y_f=array(y_c_f)
            y_f_s=smooth(y_f, int(w), window='blackman')
            y_r=array(y_c_r)
            y_r_s=smooth(y_r, int(w), window='blackman')  
            ylabel(seq1+"  (sRNA Density)")
            xlabel(ref1+"  (Reference Sequence)")
            plt.plot(x_c,y_f_s, color='r')
            plt.plot(x_c,y_r_s, color='r')
            axhline(y=0)
            plt.title('Density Map')
            #plt.ion() #Remove if running on windows
            plt.show()
        
        else:

            plt.bar(x_c,y_c_f, linewidth=1,width=1, color = 'r', edgecolor='r')
            plt.bar(x_c,y_c_r, linewidth=1,width=1, color = 'r', edgecolor='r')
            axhline(y=0, color='b')
            ylabel(seq1+"  (sRNA Density)")
            xlabel(ref1+"  (Reference Sequence)")
            plt.show()

    def den_multi_plot(self, results_list2, arg1, arg2, ref1, y_label, w=50):
        """
        Plots a density  map plot for two results lists in different colours
        
        den_multi_plot((list(str, str, int, int, str)), (list
            (str, str, int, int, str)), str, str, str, str)
        """
        a=self.org_coord()
        x_c_1= a[0]
        y_c_f_1=a[1]
        y_c_r_1=a[2]
        b=results_list2.org_coord()
        x_c_2= b[0]
        y_c_f_2=b[1]
        y_c_r_2=b[2]
        y_f_1=array(y_c_f_1)
        y_f_s_1=smooth(y_f_1, int(w), window='blackman')
        y_r_1=array(y_c_r_1)
        y_r_s_1=smooth(y_r_1, int(w), window='blackman')              
        y_f_2=array(y_c_f_2)
        y_f_s_2=smooth(y_f_2, int(w), window='blackman')
        y_r_2=array(y_c_r_2)
        y_r_s_2=smooth(y_r_2, int(w), window='blackman')  
        xlabel(ref1+"  (Reference Sequence)")
        ylabel(y_label)
        line1=plt.plot(x_c_1,y_f_s_1, color='r',label=arg1)
        line2=plt.plot(x_c_1,y_r_s_1, color='r')
        line3=plt.plot(x_c_2,y_f_s_2, color='g', label=arg2)
        line4=plt.plot(x_c_2,y_r_s_2, color='g')
        axhline(y=0)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title('Density Map')
        #plt.ion() #Remove if running on windows
        plt.show()
    

class Incomplete_Args(Exception):
    """For when not enough arguments for a particular analysis are entered into
    the GUI
    """
    
    pass  

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

def smooth(x,window_len=20,window='hamming'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]


   
run = User_Input() 
run.runner()
