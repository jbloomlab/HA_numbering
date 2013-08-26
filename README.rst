======================
HA_numbering
======================

This is a small script that converts between different numbering scheme for influenza hemagglutinin (HA). The motivation for this script is that in cloning applications it is often most convenient to number HA sequentially starting with 1 for the N-terminal methionine. However, the various PDB structures and standard naming schemes use other numbering systems. This script converts among those numbering systems.

This script was written by `Jesse Bloom`_


Installation
-----------------------

You can download the script `on GitHub`_. This is a `Python`_ script that is known to work with versions 2.6 and 2.7, and probably works with other versions 2.* as well. It has been tested on Linux and Mac OS X.

This script also requires you to install the `PROBCONS`_ alignment program, which is used to make the alignments for the number conversions. `PROBCONS`_ can be downloaded for free, and this script has been tested with version 1.12.


Numbering schemes
---------------------

This script converts between several numbering schemes. You will give it a HA protein sequence and some residue numbers in that sequence. It will then report the equivalent numbers in the following schemes:

* **sequential** is the numbering of the protein sequence that you provide in 1, 2, ... numbering. So for example, for the following sequence::

    MKAILVVLL

  the M is residue 1, the K is residue 2, etc. The HA1 and HA2 polypeptides are numbered as part of the same protein sequence in this numbering scheme.

* **4HMG** is the numbering that is used in the PDB structure `4HMG`_, which is the crystal structure of the HA from human H3N2 strain A/Aichi/2/1968 (also known as the X-31 HA). This is the numbering scheme that is often referred to as the "H3 numbering system." The HA1 and HA2 polypeptides are numbered as different sequences in this numbering scheme.

* **3TI6** is the numbering scheme that is used in the PDB structure `3TI6`_, which is the crystal structure of the HA from the human 2009 pandemic H1N1 strain A/California/4/2009. The HA1 and HA2 polypeptides are numbered as different sequences in this numbering scheme.


Input file
-------------------

To run this script, create an input file of the format below::

    probconspath /Users/jbloom/probcons/
    ha_sequence MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI
    sites 180 15 389 288 216 312 179 145

This file specifies three keys which have the following meanings:

  * *probconspath* is a path to the `PROBCONS`_ executable. Within the directory specified by this path, there should be an executable with the name ``probcons``.

  * *ha_sequence* gives the protein sequence of the HA that we are examining as a string of letters. 

  * *sites* gives the sites of interest in **sequential** numbering of the sequence specified by *ha_sequence*. This should be one or more sites indicated by integer numbers separated by spaces.



.. _`on GitHub`: https://github.com/jbloom/HA_numbering
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Python`: http://www.python.org/
.. _`PROBCONS`: http://probcons.stanford.edu/download.html
.. _`4HMG`: http://www.pdb.org/pdb/explore/jmol.do?bionumber=1&structureId=4HMG
.. _`3TI6`: http://www.pdb.org/pdb/explore/explore.do?structureId=3ti6
