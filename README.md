

# virtree

![virtree-logo](https://github.com/faylward/virtree/blob/main/virtree-logo.jpg)


virtree is a script for generating a binary alignment file of viral genomes based on the present/absence of protein families. The input is a folder of .faa files (one file for the proteins in a genome) and the result is a pseudo-FASTA file of the binary presence/absence of protein families in each genome. This pseudo-FASTA file  can then be used for building phylogenetic trees in a similar manner as a regular multi-FASTA alignment file. This was initially developed for making phylogenetic trees of bacteriophages (Caudoviricetes) where universal marker genes are lacking, but in principle this method could be applied in a range of different contexts.

### Dependencies

virtree is written in Python 3.9 and requires pandas. HMMER3 is used for hmm searches, and hmmsearch must be installed in your PATH. 

If you get weird error messages about files not existing, this would be a good first step to check. 

### Database download

virtree can be used with any HMM database - you can decide what protein family HMMs make the most sense for your analysis. 

A pre-compiled set of VOGs that we use for treebuilding of bacteriophages can be downloaded and set up using the code below:

> wget -O hmm.tar.gz  <https://zenodo.org/records/12744518/files/hmm.tar.gz?download=1>

and then

> tar -xvzf hmm.tar.gz

This should create a hmm/ directory with the vog_05p.hmm and vog_025p.hmm files that can be used for analyzing Caudoviricetes. 
The only difference between these databases is that the vog_025p.hmm file has a few more models and takes a bit longer to run. 

### Basic Usage

To test if hmmsearch will run properly type: \

> python virtree.py -i test_proteins -p test_profile -t 2

This should take ~30 seconds to run and create the following files:

\*_bin.fasta - this is a binary FASTA file - i.e. a FASTA formatted file with 0's or 1's for the absence or presence of certain protein families. The order of the protein families is the same for each record. 

\*._profile.tsv - this is a tab-delimited file with the presence/absence of the protein families. You can load this into R if you wish to check the occurrence of specific protein families. 


### Options

**-db, --database** This is the PATH to the HMM database that is used. This could be any database you wish - Pfam, COG, etc. Just make sure to provide the full PATH. The default is hmm/vog_05p.hmm, so it is assumed that you are running this while located in the virtree directory. 

**-g, --minhit** Minimum number of hits against the HMM database that must be recorded if a genome is to be included in the output binary fasta. Default is 2. 

**-e, --evalue** E-value to use for hmmsearch - default is 1e-3

**-t, --cpus** how many CPUs to use for hmmsearch - default is 1


                                                                                                                        

                                  

                                                                 
                               *+*                               
                           %=+-----=-                            
                        *===---------=--*                        
                      %=#+#===========*-+.#                      
                      @=+#-------------+:.                       
                      @==-=-----------=::.                       
                      @==-------------::-.                       
                      @=-------------::::.                       
                      @=----+--------::::.                       
                      @*-----=-----=::::::                       
                      @------------::::::=                       
                      @-------=---::::::::                       
                      @#-------=-=::::::-#%                      
                          *==+===--:..=                          
                             %===..                              
                               #*#                               
                            ==*+=-++%                            
                           = -+++=+-*                            
                          %# *=---=+ =                           
                          =  -*++=+- #%                          
                         +  ====+==+% =                          
                        %*  =------== @%                         
                            -+=---+--                            
                             +=+*+-+                             
                             -=====-                             
                             =-----=                             
                      ++   --*----=:     *=                      
                +#   *  = %-%+=*+=-=    + #@   *-                
               =  %= =   +#= -==+==-   *   + =@  :               
              =     @=@   *+ *+----- @*   @=@     =%             
             *      +  =@ ** #=++==+ *  %+  *      *@            
           %#      @%    +=%===*+-=-% %*     %      ##           
          %%       *      -**--++=--*#       *       %*          
         +         %      -  **=+=++          %        -         
        =                 -*+*+---+--+                  :        
       -                 @-==++++++:::=                  :       
                         #+*=------=:*                           
                         +#                                      
                         =@                                      
                         -                                       
                                                                 
                                                                 

                                                      
                                                                                                                                            

