#!perl -w

# this program removes those absolutely fucking annoying 
# ^M characters at the end of the line in files generated
# by MSword and Excel
# usage: perl remcarM.pl filename
# caution - this erases the original version of the file


$file = shift(@ARGV) or die "error: no infile\n";
system(qq{perl -p -i -e 's/\015/\n/go' $file});
system(qq{perl -p -i -e 's/\015\012/\n/go' $file});
