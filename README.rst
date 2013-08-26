======================
HA_numbering
======================

Summary
-----------

This is a small script (*HA_numbering.py*) that converts between different numbering scheme for influenza hemagglutinin (HA). The motivation for this script is that in cloning applications it is often most convenient to number HA sequentially starting with 1 at the first residue. However, various PDB structures and publications schemes use other numbering systems. This script converts among some of those numbering systems.

This script was written by `Jesse Bloom`_


Installation
-----------------------

You can download the script `on GitHub`_ (just click on the button that says ``Download ZIP`` on the right side of the page partway down). This is a `Python`_ script that is known to work with versions 2.6 and 2.7, and probably works with other versions 2.* as well. So you will need to have `Python`_ installed on your computer. 

This script also requires you to install the `PROBCONS`_ alignment program, which is used to make the alignments for the number conversions. `PROBCONS`_ can be downloaded for free, and this script has been tested with version 1.12.


Numbering schemes
---------------------

This script converts between several numbering schemes. You will give it a HA protein sequence and some residue numbers in **sequential** numbering of that sequence. It will then report the equivalent numbers in the following schemes:

* **sequential** is the numbering of the protein sequence that you provide in 1, 2, ... numbering. So for example, for the following sequence::

    MKAILVVLL

  the M is residue 1, the K is residue 2, etc. The HA1 and HA2 polypeptides are numbered as part of the same protein sequence in this numbering scheme.

* **4HMG** is the numbering that is used in the PDB structure `4HMG`_, which is the crystal structure of the HA from human H3N2 strain A/Aichi/2/1968 (also known as the X-31 HA). This is the numbering scheme that is often referred to as the "H3 numbering system." The HA1 and HA2 polypeptides are numbered as different sequences in this numbering scheme.

* **4JTV** is the numbering scheme that is used in the PDB structure `4JTV`_, which is the crystal structure of the HA from the human 2009 pandemic H1N1 strain A/California/4/2009. The HA1 and HA2 polypeptides are numbered as different sequences in this numbering scheme.


Input file
-------------------

To run this script, create an input file of the format below::

    probconspath /Users/jbloom/probcons/
    ha_sequence MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI
    sites 180 15 389 288 216 312 179 145

This file specifies three keys which have the following meanings:

  * *probconspath* is a path to the `PROBCONS`_ executable. Within the directory specified by this path, there should be an executable with the name ``probcons``. Right now you must specify a valid directory here even if `PROBCONS`_ is in the current search path.

  * *ha_sequence* gives the protein sequence of the HA that we are examining as a string of letters. 

  * *sites* gives the sites of interest in **sequential** numbering of the sequence specified by *ha_sequence*. This should be one or more sites indicated by integer numbers separated by spaces.


Running the script
-------------------

To run the script, create an input file of the format described above and put it into the same directory as the *HA_numbering.py* script, and then run the script with the input file as the sole argument. For example, if you name your input file *example_infile.txt* (you can name it anything you want), you would then run::

    python HA_number.py example_infile.txt

If you formatted the input file correctly, specified valid numbers for the *sites* variable in your input file for the HA sequence that you provided for the *ha_sequence* variable, and provided a valid value for *probconspath*, then the program should print output giving the site number mappings.


Output format
---------------

The output format is printed to standard output. For example, for the example input file given above, you should get the following output::

    Beginning execution of HA_numbering.py script.
    Reading input from example_infile.txt
    Making PROBCONS alignments...
    Alignments complete.

    Here are the corresponding residue numbers:

    site K180 in sequential numbering of your HA sequence corresponds to:
      * V166 in HA1 in 4HMG
      * K169 of HA1 in 4JTV

    site A15 in sequential numbering of your HA sequence corresponds to:
      * N8 in HA1 in 4HMG
      * an alignment gap in 4JTV

    site I389 in sequential numbering of your HA sequence corresponds to:
      * I45 in HA2 in 4HMG
      * I45 of HA2 in 4JTV

    site P288 in sequential numbering of your HA sequence corresponds to:
      * P273 in HA1 in 4HMG
      * P277 of HA1 in 4JTV

    site V216 in sequential numbering of your HA sequence corresponds to:
      * V202 in HA1 in 4HMG
      * V205 of HA1 in 4JTV

    site I312 in sequential numbering of your HA sequence corresponds to:
      * V297 in HA1 in 4HMG
      * I301 of HA1 in 4JTV

    site S179 in sequential numbering of your HA sequence corresponds to:
      * N165 in HA1 in 4HMG
      * S168 of HA1 in 4JTV

    site S145 in sequential numbering of your HA sequence corresponds to:
      * Q132 in HA1 in 4HMG
      * S134 of HA1 in 4JTV

    Script complete.



.. _`on GitHub`: https://github.com/jbloom/HA_numbering
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Python`: http://www.python.org/
.. _`PROBCONS`: http://probcons.stanford.edu/download.html
.. _`4HMG`: http://www.pdb.org/pdb/explore/jmol.do?bionumber=1&structureId=4HMG
.. _`4JTV`: http://www.pdb.org/pdb/explore/explore.do?structureId=4jtv
