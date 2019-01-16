# BT_NF-kB
The Scripts used in my bachelor thesis "Bioinformatic Prediction of NF-kB Target Genes in the Human Genome."
Resulting datasets are not put online, as they are intended for internal use.

For a guide on how to install third party python modules (listed below), which the scripts need to work properly, you can check https://automatetheboringstuff.com/appendixa/ . This ressource is taken from the book "Automate the boring stuff with python" (available for free on https://automatetheboringstuff.com/#toc ), by Al Sweigart, my main ressource for learning the Python used in this thesis.  

Dependencies:

:: OS: Everything should run on Windows and MacOS(tested only partially), Linux not tested.

:: Firefox must be installed

:: Third party modules: 

                        - re

                        - requests
                        
                        - csv
                        
                        - sys
                        
                        - os
                        
                        - selenium
                        
                        - time
                        
                        - lxml
                        
                        - platform
                        
                        - numpy
                        
                        - matplotlib
                        
                        - scipy
                        
                        - itertools
                        
                        
:: Two Folders in the scripts directory (names mandatory):

    - "geckodriver" -> holds the geckodriver utility for your OS (downloadable at https://github.com/mozilla/geckodriver/releases )
    
    - "source_files" -> must hold: 
    
                                   - "gene2ensembl" mapping file from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
    
                                   - "cross_references.txt" & "promoter_ensembl.txt" from ftp://ccg.vital-it.ch/epdnew/H_sapiens/006/db/
                                   
                                   - "all_data.txt" and "promoters.gff3" the result of the scripts by bpucker on https://github.com/bpucker/NFkB
                                   - "fasta_sequences_from_-5000_to_+100.fa" downloadable via EPD's select/download tool at https://epd.vital-it.ch/EPDnew_select.php by not entering any ids, pressing select and specifying a range of -5000 to +100 on the website after. Name needs to be adjusted after downloading.
                                   
::"Figure_6.py" MUST be executed FIRST and requires some manual steps, which are described in the scripts comments and printed to the console.

:: Most of the scripts depend on the Bielefeld_basic_python_dicts.py" script, so this must be downloaded and put in the directory the script you want to execute is in.
