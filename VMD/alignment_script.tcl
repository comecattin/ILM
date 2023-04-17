

# Load trajectory and reference structure
mol new /home/ccattin/Documents/Cluster/278K/ES/clustering/C1-f12652-s4327.pdb type pdb

# Create new molecule representation for protein
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 2
mol addrep top

lappend auto_path /home/ccattin/Documents/Code/VMD/la1.0
lappend auto_path /home/ccattin/Documents/Code/VMD/orient

package require Orient
namespace import Orient::orient

set sel [atomselect top "all"]
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $sel]



proc rotate_axis {vec deg {molid top}} {
    # get the current matrix
    lassign [molinfo $molid get rotate_matrix] curr
    # the transformation matrix
    set r [trans axis $vec $deg]
    # get the new matrix
    set m [transmult $r $curr]
    # and apply it to the molecule
    molinfo $molid set rotate_matrix "{ $m }"
}

graphics top delete all

#rotate_axis {0 1 0} 90