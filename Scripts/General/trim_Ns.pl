use warnings;
use strict;
my %hashR;


open(IN, "/bin/gzip -dc $ARGV[0] |");
#print $ARGV[0];
#open(DEB, ">$ARGV[]$ARGV[0].debug.txt");
$ARGV[0] =~ s/.fastq.gz//;
open (DEB, "| /bin/gzip -c > $ARGV[1]/$ARGV[0].debug.gz") or die "error starting gzip $!";
open (FIVEP, "| /bin/gzip -c > $ARGV[1]/$ARGV[0].5p.gz") or die "error starting gzip $!";
open (OUT, "| /bin/gzip -c > $ARGV[1]/$ARGV[0].trimN3p.fastq.gz") or die "error starting gzip $!";
#open(OUT, ">$ARGV[0].txt");

my $c = 0;
while (my $line = <IN>){
  $c++;
  #print "$c boia \n";
  chomp $line;
  $hashR{$c} = $line;
  if ($c == 4){
    my $lenN = 0;
    #print "$hashR{2} ccio \n";
    if ($hashR{2} =~ /(N+)$/){
      $lenN = length($1);
      #print "$lenN lenn \n";
      my @split2 = split('',$hashR{2});
      splice(@split2, scalar(@split2)-$lenN, $lenN);
      $hashR{2} = join("", @split2);

      my @split4 = split('',$hashR{4});
      splice(@split4, scalar(@split4)-$lenN, $lenN);
      $hashR{4} = join("", @split4);
    }

    if ($hashR{2} =~ /^(N+)/){
        my $lenN5p = length($1);
        print FIVEP "$hashR{1}\t$lenN5p\n";  
    }

    my $g =0;
    while ($g < 4){
      $g++;
      print OUT "$hashR{$g}\n";
    }
    print DEB "$hashR{1}\t$lenN\n";
    
    $c=0;
  }
}

close(IN);
close(OUT);
close(DEB);
close(FIVEP);

exit
