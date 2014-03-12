#!/bin/sh
# Simon Warfield simon.warfield@childrens.harvard.edu
# $Id:$
# \
exec tclsh $0 ${1+"$@"}

if {($argc != 8)} {
  puts "renumberFiles.tcl inPrefix firstfile lastfile inc "
  puts " outPrefix outFirst outLast outInc"
  puts "Example: renumberFiles.tcl I.%03d 40 50 1 new.%03d 30 40 1"
  puts "The argument count argc is $argc"
  exit 1
}

set inPrefix [lindex $argv 0]
set infirst [lindex $argv 1]
set lastfirst [lindex $argv 2]
set inc [lindex $argv 3]

set outprefix [lindex $argv 4]
set outfirst [lindex $argv 5]
set outlast [lindex $argv 6]
set outinc [lindex $argv 7]

for {set i $infirst} {$i <= $lastfirst} {incr i $inc} {
  set inname [format ${inPrefix} $i]
  set outnumber [expr $outfirst + ($i - $infirst)*$outinc/$inc]
  set outname [format ${outprefix} $outnumber]
puts "inname $inname i $i outnumber $outnumber outname $outname"
  exec cp [format ${inPrefix} $i] ${outname}
}

exit 0
