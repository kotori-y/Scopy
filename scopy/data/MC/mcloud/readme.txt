
Molecule Cloud - Instructions

Compile the program using the command:
javac ertl/mcloud/MCloud.java

Run it with:
java ertl/mcloud/MCloud -f scaffolds.smi

This default version displays only rectangles instead of molecular structures. To add chemical intelligence, implement the MolHandler interface with your favorite cheminformatics toolkit for parsing SMILES and rendering molecules.

A MolHandler implementation using the free open source Avalon Cheminformatics Toolkit is included in this distribution (AvalonMolHandler.java).

Below are instructions for building and running Molecule Cloud with the Avalon Toolkit on a Windows 32-bit operating system:

1 - modify the source of MCloud.java (around line 570), comment the call to the constructor of TestMolHandler and uncomment the one for AvalonMolHandler 

2 - download Avalon Cheminformatics Toolkit from SourceForge
http://sourceforge.net/projects/avalontoolkit/

3 - copy the following files from the Avalon Toolkit distribution (they are located in the directory libraries/) into the mcloud directory, i.e. the directory from which you want to run the program:
depictjni.jar
depict.jar
avalon_jni_JNISmi2Mol.dll

4 - compile as:
javac -cp ".;depictjni.jar;depict.jar" ertl/mcloud/MCloud.java ertl/mcloud/AvalonMolHandler.java

5 - run as:
java -cp ".;depictjni.jar;depict.jar" ertl/mcloud/MCloud -f scaffolds.smi -n 150
(using the most common scaffolds from ZINC as data)

To use the program on another operating systems consult the Avalon Toolkit documentation.
For any questions concerning the Avalon Toolkit contact its author Bernhard Rohde.

To see Molecule Cloud options, run the program without any options, or check the source code.

Good luck!


Molecule Cloud theory is available in:
P. Ertl, B. Rohde: The Molecule Cloud - compact visualization of large collections of molecules
J. Cheminformatics 2012, 4:12
http://www.jcheminf.com/content/4/1/12/

If you are using the code, cite please the article.
