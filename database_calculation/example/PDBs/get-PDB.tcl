set n [molinfo top get numframes]

set all [atomselect top "protein or segname MEMB"]
for {set i 0} {$i < $n} {incr i} {
    $all frame $i
    $all writepdb frame_$i.pdb
}
