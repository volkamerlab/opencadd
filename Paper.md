---
title: 'StructuralAlignment: A Python library to use different methods for structural alignment.'
tags:
    - python
    - bioinformatics
authors:
 - name: Julian Pipart
   orcid: 0000-0002-1089-5931
 - name: Jaime Rodríguez-Guerra
 - name: Andrea Volkamer
 - name: Dennis Köser
 - name: Annie Pham
 - name: Enes Kuni
date: 29 April 2020
---
#StructuralAlignment
##A Python library to use different methods for structural alignment.

###Authors
- Julian Pipart
- Jaime Rodríguez-Guerra
- Andrea Volkamer
- Dennis Köser
- Annie Pham
- Enes Kuni

# Abstract

**Motivation**: There are a lot of different methods to align structures with each other. But most of the time you do not know which method is the best for a specific situation, or you want to build a python script around these methods. Some methods are not even implemented in python. That makes it difficult to import these methods in custom scripts or software projects. In order to use these different methods you need to install them separately. They also return output which is not consistent over the different methods.
**Results**: StructuralAlignment has been developed to bring multiple methods for structural alignment together and provide consistency in the output as well as input for any Python 3.7 or Python 3.8 interpreter. It provides the possibility of interactive programming and makes it easier to use third-party software for structural alignment.
**Availability and implementation**: StructuralAlignment is MIT-licensed and available at https://github.com/volkamerlab/StructuralAlignment
**Contact**: julian.pipart@fu-berlin.de or jaime.rodriguez@charite.de


# 1 Introduction
Python is a very powerful programming language. No matter if it is used for projects in big companies or for small scripts. It is especially useful for solving scientific or computational problems because of its variety of different packages to use. This includes problems of structural bioinformatics and chemistry. There is already a lot of software to solve specific problems and most of them use the same idea but a different approach. They often accept not the same input and provide different outputs. This can be really hindering for the research or the process of software development.


# 2 Materials and methods
StructuralAlignment is a Python package which brings multiple methods for structural alignment together. The user can specify the method to use and does not have to worry about the dependencies and different input parameters for different methods. To ensure consistency across all methods StructuralAlignment uses Atomium ([Ireland, S. M., & Martin, A. C., 2020](#references)). Atomium is able to catch the structures the user gives as input from the RCSB. With the help of Atomium StructuralAlignment first transforms the structure to an atomium.model and passes this model on to the method of choice. All methods calculate a transformation matrix with their own approaches. This matrix is applied to the original model. By that there will be no metadata lost. The output of StructuralAlignment is also an atomium.model for all methods implemented. The model can be saved as a Protein Data Bank file (PDB) or similar filetypes. Because most methods do not know what an atomium.model is, they save the model as a temporary PDB file and parse it. This makes the package very easy to use in development of software. It also assures the ability to be further extended. The correctness is ensured by the Github CI workflow. If a method is written in a different programming language and is open source, it is wrapped in Python. Other methods, like Chimera matchmaker ([Pettersen, E. F., Goddard, T. D., Huang, C. C., Couch, G. S., Greenblatt, D. M., Meng, E. C., & Ferrin, T. E. 2004](#references)) has been reimplemented in this library.

# 3 Results
StructuralAlignment is a command line application with support for the latest Linux distributions (it is tested for Python 3.7 and Python 3.8 on Ubuntu 18.04) and Mac OS X (it is tested for Python 3.7 and Python 3.8 on Mac OS X 10.15.4). By installing StructuralAlignment, all dependencies needed are installed with it. Only Python 3.7 or later is required before installation. The installation is made possible through conda with the command: `conda install StructuralAlignment`. It is the same command for both platforms. A great advantage is, that StructuralAlignment uses Python 3. That guarantees good compatibility with most modern software written in Python.

**3.1 Usage**
Running StructuralAlignment in the terminal can be performed with the command: `structuralalignment TARGET_STRUCTURE STRUCTURE_TO_ALIGN --method METHOD_OF_CHOICE --method-options METHOD_OPTIONS`. StructuralAlignment features the methods MMLigner ([Collier, J. H., Allison, L., Lesk, A. M., Stuckey, P. J., Garcia de la Banda, M., & Konagurthu, A. S., 2017](#references)), Theseus ([Theobald, D. L., & Wuttke, D. S., 2006](#references)) and Chimera's MatchMaker. It is also possible to perform multiple pairwise alignments by adding more structures to align to the command above. All structures will therefore be aligned against the target structure. The output in the terminal is the RMSD.

**3.2 Demonstrating the Benefits of 'StructuralAlignment'**
StructuralAlignment brings different methods for structural alignment together, so that it is not necessary to have all the stand alone programs installed. The methods also have different approaches.
For the more traditional approach, one can use the method matchmaker from chimera
```
INPUT: structuralalignment 2GZ9 5R8T --method=matchmaker
OUTPUT:
  _____ _                   _                   _
  / ____| |                 | |                 | |
 | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ __ _| |
  \___ \| __| '__| | | |/ __| __| | | | '__/ _` | |
  ____) | |_| |  | |_| | (__| |_| |_| | | | (_| | |
 |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_|
     /\   | (_)                                | |
    /  \  | |_  __ _ _ __  _ __ ___   ___ _ __ | |_
   / /\ \ | | |/ _` | '_ \| '_ ` _ \ / _ \ '_ \| __|
  / ____ \| | | (_| | | | | | | | | |  __/ | | | |_
 /_/    \_\_|_|\__, |_| |_|_| |_| |_|\___|_| |_|\__|
                __/ |
               |___/

 v0+untagged.91.g31533ad · Brought to you by @volkamerlab

RMSD for alignment between `2GZ9` and `5R8T` is 0.4Å
```

If one prefer a statistical approach, one can use for example Theseus.
```
INPUT: structuralalignment 2GZ9 5R8T --method=theseus
OUTPUT:
  _____ _                   _                   _
  / ____| |                 | |                 | |
 | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ __ _| |
  \___ \| __| '__| | | |/ __| __| | | | '__/ _` | |
  ____) | |_| |  | |_| | (__| |_| |_| | | | (_| | |
 |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_|
     /\   | (_)                                | |
    /  \  | |_  __ _ _ __  _ __ ___   ___ _ __ | |_
   / /\ \ | | |/ _` | '_ \| '_ ` _ \ / _ \ '_ \| __|
  / ____ \| | | (_| | | | | | | | | |  __/ | | | |_
 /_/    \_\_|_|\__, |_| |_|_| |_| |_|\___|_| |_|\__|
                __/ |
               |___/

 v0+untagged.91.g31533ad · Brought to you by @volkamerlab

RMSD for alignment between `2GZ9` and `5R8T` is 1.7Å
```


**3.3 Current limitations**
For now, the only method, that supports a benchmark is MMLigner in form of the ivalue. The other methods do not provide a benchmark yet. It is not clear if the ivalue is a good and consistent benchmark over different methods. This means it is not obvious whether an alignment is good or not, compared to the other alignments. It is also not possible to perform a multiple alignment at the moment. One can perform multiple pairwise alignments successively.

**3.4 Future work**
In the future more alignment methods will be added. Because StructuralAlignment is open source, every developer can participate in making StructuralAlignment better and more versatile. The MMLigner method already provides a benchmark with the ivalue. The developer team is working on providing a benchmark for every alignment method. This benchmark should be consistent over the different methods. In addition to the benchmarking the team is working on a solution to visualize the aligned structures in this library, without the need of any third-party program.

# Acknowledgments
The authors want to thank the team of Volkamerlab, and more particularly Jaime Rodríguez-Guerra and Andrea Volkamer for making this project possible and for leading the development process.

# References
Collier, J. H., Allison, L., Lesk, A. M., Stuckey, P. J., Garcia de la Banda, M., & Konagurthu, A. S. (2017). Statistical inference of protein structural alignments using information and compression. Bioinformatics, 33(7), 1005-1013.

Ireland, S. M., & Martin, A. C. (2020). atomium—a Python structure parser. Bioinformatics.

Theobald, D. L., & Wuttke, D. S. (2006). THESEUS: maximum likelihood superpositioning and analysis of macromolecular structures. Bioinformatics, 22(17), 2171-2172.

Pettersen, E. F., Goddard, T. D., Huang, C. C., Couch, G. S., Greenblatt, D. M., Meng, E. C., & Ferrin, T. E. (2004). UCSF Chimera—a visualization system for exploratory research and analysis. Journal of computational chemistry, 25(13), 1605-1612.

