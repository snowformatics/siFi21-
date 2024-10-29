siFi: Software for long double-stranded RNAi-target design and off-target prediction

Version v1.2.3-0008

Software Download: [http://www.snowformatics.com/si-fi.html](https://sourceforge.net/projects/sifi21/files/sifi1.2.3-0008.exe/download)

Reference: 
URL=https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2019.01023

DOI=10.3389/fpls.2019.01023

Stefanie Lück, Tino Kreszies, Marc Strickert, Patrick Schweizer, and Dimitar Douchkov

RNA interference (RNAi) is a nucleic acid complementarity-based biological phenomenon and a widespread natural mechanism for the inhibition or regulation of gene expression. RNAi is an essential part of the immune response to viruses and other foreign genetic material especially in plants but also in many fungal and animal species where it is a part of their innate immunity system. RNAi has become an important research tool for studying gene function by strong and selective suppression of the genes of interest also in large-scale screens. However, the application of RNAi as a technology raises important questions about efficiency and specificity of corresponding gene constructs. Since the RNAi machinery selects targets based on sequence similarity, there is an inherent risk of hitting non-targeted genes. Therefore a reliable search and prediction of off-targets is crucial for the practical application of RNAi. Besides being gene-specific, a successful RNAi construct should also be able to produce a strong silencing effect, which is highly depending on selecting an optimal sequence for design of the RNAi construct.
Here, we report on a software called siFi for optimizing long double-stranded RNAi- target design and for prediction of RNAi off-targets. It is open source desktop software that provides an intuitive graphical user interface, works in Microsoft Windows environment and can use custom sequence databases in  standard FASTA format. 

Attribution-NonCommercial-ShareAlike 2.0 Generic (CC BY-NC-SA 2.0) License
# Installation and requirements
### requirements
- Anaconda
- bowtie
- FASTA file to check against (Ensembl cDNA for example)

### Installation
```
conda create --name sifi python=3.9
conda activate sifi
conda install --yes --file requirements.txt
```

### Running
`python main.py`

# We have received many requests to save the siRNA sequences, and we understand the importance of this feature
Unfortunately, this function has not been implemented, and si-Fi is currently not updated due to the need for a complete rewrite to meet current Python standards and library versions. 

However, there is a workaround, though it requires a bit of effort. You will need to locate the system's temporary directory, which on Windows is typically:

C:\Users\YOUR_NAME_HERE\AppData\Local\Temp

After running si-Fi, you will find a temporary file with a .json extension, such as tmppaivxu.json.

To locate this file:

- Check the time of the run.
- Sort the files in the temporary directory by the latest creation time.
  
This file contains all the siRNAs along with the relevant information in JSON format. There are numerous JSON parsers available online that can help you convert this data into a table-like format.

We apologize for the inconvenience, but currently, this is the only method available to extract the siRNA sequences.
