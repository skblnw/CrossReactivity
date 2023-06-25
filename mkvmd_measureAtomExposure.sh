#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL_PEPTIDE="segname PROC"
SEL_pHLA="segname PROA PROC"
SEL_HLA_HELIX="segname PROA and resid 56 to 85 138 to 174"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2"; exit 1; }

files=("$PDB" "$TRJ")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

 # and not name CA C N O HN HA

echo "O,N,CS,H," > ${OUTPUT}_area_peptide
echo "O,N,CS,H," > ${OUTPUT}_area_hlahelix
cat > tcl << EOF
proc measureAtomExposure1 { nn } {
    set total_area_buried 0
    set output []

    set sel [atomselect top "$SEL_PEPTIDE" frame \$nn]
    set selall [atomselect top "$SEL_pHLA" frame \$nn]

    foreach element {{$SEL_PEPTIDE and name O "O.*"} {$SEL_PEPTIDE and name N "N.*"} {$SEL_PEPTIDE and name C S "C.*" "S.*"} {$SEL_PEPTIDE and name H "H.*"}} {
        set rest [atomselect top \$element frame \$nn]
        set area_exposed [format "%.2f" [expr [measure sasa 1.4 \$selall -restrict \$rest]]]
        lappend output \$area_exposed
        puts -nonewline "[format "%3.0f" \$area_exposed] + "
        incr total_area_buried [format "%.0f" \$area_exposed]
    }

    puts "Total Exposed=\$total_area_buried"
    puts "Total Residue=[format "%3.0f" [measure sasa 1.4 \$selall -restrict \$sel]]"
    return \$output
}

proc measureAtomExposure2 { nn } {
    set total_area_buried 0
    set output []

    set sel [atomselect top "$SEL_HLA_HELIX" frame \$nn]
    set selall [atomselect top "$SEL_pHLA" frame \$nn]

    foreach element {{$SEL_HLA_HELIX and name O "O.*"} {$SEL_HLA_HELIX and name N "N.*"} {$SEL_HLA_HELIX and name C S "C.*" "S.*"} {$SEL_HLA_HELIX and name H "H.*"}} {
        set rest [atomselect top \$element frame \$nn]
        set area_exposed [format "%.2f" [expr [measure sasa 1.4 \$selall -restrict \$rest]]]
        lappend output \$area_exposed
        puts -nonewline "[format "%3.0f" \$area_exposed] + "
        incr total_area_buried [format "%.0f" \$area_exposed]
    }

    puts "Total Exposed=\$total_area_buried"
    puts "Total Residue=[format "%3.0f" [measure sasa 1.4 \$selall -restrict \$sel]]"
    return \$output
}

# /------------------/
# /     Main Body    /
# /------------------/

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      set outf1 [open ${OUTPUT}_area_peptide "a"]
      set outf2 [open ${OUTPUT}_area_hlahelix "a"]

      # Call calc funtion you like
      set output1 [measureAtomExposure1 \$nn]
      set output2 [measureAtomExposure2 \$nn]

      # Write to file
      foreach element \$output1 {puts -nonewline \$outf1 "\$element,"}
      puts \$outf1 ""
      foreach element \$output2 {puts -nonewline \$outf2 "\$element,"}
      puts \$outf2 ""
      # Remember to close the file
      close \$outf1
      close \$outf2
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
