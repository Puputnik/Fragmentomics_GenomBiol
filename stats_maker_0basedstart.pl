use warnings;
use Storable;
use List::Util qw[min max];

my %hash;

#### parsing input sam
while (<STDIN>){
    my $line = $_;
    
    unless ($line =~ /^@/){
        my @split = split(/\t/, $line);
        
        unless ( $split[2] eq '*'){
        
        #### calculating read length from cigar (meaningful only for Nanopore)
        my @cig = ( $split[5] =~ /(\d+\w)/g);
        my $len = 0;
        for my $c (@cig){
            if ($c =~ /(\d+)[MX=DN]/g){
                $len = $len + $1;
            }
        }
        
        my $cigo = $split[5];
        
        my $H1 = 0;
        my $H2 = 0;
            
        my $S1 = 0;
        my $S2 = 0;
        
        
        #### calculating length of leftmost hardclip
        if ($cigo =~ m/^(\d+)H/){
            $H1 = $1;
            $cigo =~ s/^$H1\H//;
        } 
        
        #### calculating length of rightmost hardclip
        if ($cigo =~ m/(\d+)H$/){
            $H2 = $1;
            $cigo =~ s/$H2\H$//;
        }
        
        #### calculating length of leftmost softclip
        if ($cigo =~ m/^(\d+)S/){
            $S1 = $1;
        } 
        
        #### calculating length of rightmost softclip
        if ($cigo =~ m/(\d+)S$/){
            $S2 = $1;
        
        }
     
        my $readid  =  $split[0];
        my $flag    =  $split[1];
        my $contig  =  $split[2];
        my $start   =  $split[3] - 1 ; #### converting read start from 1-based to 0-based (for downstream getfasta tool)
        my $end     =  $start + $len ;
        my $Q       =  $split[4];
        
        #### getting mapping strand from flag
        my $strand;
        if ($flag & 0x10){
            $strand  =  "-";
        } else {
            $strand  = "+";
        };
        
        #### getting mate information (forward or reverse) from flag (meaningful only for Illumina)
        my $mate;
        if ($flag & 0x40){
            $mate  =  "F";
        } elsif ($flag & 0x80){
            $mate  = "R";
        } else {
            $mate = "NA";
        };
        
        #### getting TLEN field (meaningful only for Illumina)
        my $TLEN       =  $split[8];
        
        
        #### printing
        print $contig ; print "\t" ;
        print $start  ; print "\t" ;
        print $end    ; print "\t" ;
        print $readid ; print "\t" ;
        print "0\t";
        print $strand ; print "\t" ;
        print $flag   ; print "\t" ;
        print $Q      ; print "\t" ;
        print $H1     ; print "\t" ;
        print $S1     ; print "\t" ;
        print $len    ; print "\t" ;
        print $S2     ; print "\t" ;
        print $H2     ; print "\t" ;
        print $mate   ; print "\t" ;
        print $TLEN   ; print "\n" ;
        
        }   
    }
}

exit
