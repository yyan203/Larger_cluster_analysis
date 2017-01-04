#!/usr/bin/perl
use Shell;

#######################################################################
# A Perl script example of importing data and conduct simple analysis #
#######################################################################


####### Command line parameter processing####################################

if($#ARGV < 6) {
  print "Correct syntax: lmps.analysis.pl inputfile xdif atomtype1 atomtype2 ghosttype starttimestep lowestY(to reduce the calculation burden) (default frequency:5 10 20 50,modify code:frequencylist, to suit your need)\n";
  print "Attention!!! atom type is substracted by 1 to suit cluster analysis using MapCorrelation3 28\n";
  exit;
}

#my ($inputfile, $outputfile, $inputparameter) = @ARGV;
$inputfile = shift(@ARGV);
$inputparameter = shift(@ARGV);
$inputparameter2 = shift(@ARGV);
$inputparameter3 = shift(@ARGV);
$inputparameter4 = shift(@ARGV);
$inputparameter5 = shift(@ARGV);
$inputparameter6 = shift(@ARGV);

#if($ARGV>=3) {$nextpara = shift(@ARGV);}
#else {$nextpara = 1;}

print "$inputfile\n";
print "Distance Cut:$inputparameter\n";

#my @frequencylist=(1,5,10,12,20,50);
#my @frequencylist=(12);
my @frequencylist=(1);
#my @frequencylist=(5,10,20,50);


foreach (@frequencylist) {
my $frequency=$_;

##### Import the data file#################################################
my $typecount;
my @x,@y,@z1,@z2;
my @check;
my $iframe=0;
my $i,$j;
my $flag1=0;
my $ndebris=0;
my $natoms=100000;
my $line;
open(DATA, "$inputfile");




while(<DATA>) {
   $typecount=0;


  chomp;
  if(/ITEM: TIMESTEP/) { # if($_ =~ m/ITEM: TIMESTEP/)
    $line = <DATA>; chomp($line); $timestep=$line; # actual timesteps
    $line = <DATA>; chomp($line);                  #ITEM: NUMBER OF ATOMS
    $line = <DATA>; chomp($line); $natoms=$line;   #number of atoms
    $line = <DATA>; chomp($line);                  # ITEM: BOX BOUND
    $line = <DATA>; chomp($line);                  # xlo xhigh
                                  ($xlo, $xhi)=split(/\s+/, $line);
    $line = <DATA>; chomp($line);                  # ylo yhigh
                                  ($ylo, $yhi)=split(/\s+/, $line);
    $line = <DATA>; chomp($line);                  # zlo zhigh
                                  ($zlo, $zhi)=split(/\s+/, $line);
    $line = <DATA>; chomp($line);                  # ITEM: ATOMS
    #print "$natoms $timestep\n";
    #print "$xlo, $xhi, $ylo, $yhi, $zlo, $zhi\n";
    

    if( $timestep>=$inputparameter5 && ($timestep-$inputparameter5)==int(($timestep-$inputparameter5)/($frequency*1000))*$frequency*1000 ) {
	print "Timestep is: $timestep\n";
    my $k=0; 
    my $lx=0; 
    my $lz=0; 
    my $m=0;
    my $n=0;
    
    if($timestep>$inputparameter5){open(OLIVIA, ">olivia.cutoff.$inputparameter.frequency.$frequency.$iframe");
	$Lx=$xhi-$xlo;$Ly=$yhi-$ylo;$Lz=$zhi-$zlo;
	print OLIVIA "$Lx $Ly $Lz 0.000000\n";
	}	

     for $i (0..$natoms-1) {
      #loop the information of every atom
      $line = <DATA>; chomp($line);
      ($id, $itype, $ix, $iy, $iz) = split(/\s+/, $line);
      #print ("$ix, $iy, $iz\n");
      #Store the information into arrays for later use
      # store atom type into x[][][3]

	if($timestep==$inputparameter5)
	{
	   $debriflag[$id]=0;
           $refx[$id][0] = $ix;
	   $refx[$id][1] = $iy;
	   $refx[$id][2] = $iz;
	}

	elsif ($timestep>$inputparameter5) {
	
	$x2 = $ix-$refx[$id][0];
	$y2 = $iy-$refx[$id][1];
	$z2 = $iz-$refx[$id][2];
	$displacement = sqrt($x2*$x2+$y2*$y2+$z2*$z2);

	if($debriflag[$id]==0 && ($displacement>$inputparameter) && ($itype<3) && ($iy>$inputparameter6)) {
	
	if($id==117541) {print "117541 atom displacement:$displacement\n";}

	$debriflag[$id]=1;
	$itype -= 1;
	$ndebris++;
	print OLIVIA "$id $itype $ix $iy $iz\n";} 
	}	       
      
	} ## for loop atoms

     ## delete olivia file if it is empty. etc. no debris atom generated.
     if($timestep>$inputparameter5){close(OLIVIA);
     open(OLIVIA, "olivia.cutoff.$inputparameter.frequency.$frequency.$iframe");
     $line=<OLIVIA>;print "$line"; ## first line 
     $line=<OLIVIA>;print "$line"; ## first line 
     if($line){
      		close(OLIVIA);
      		 }
     else {
     	    print "olivia file is empty, there is no debris atom in this period, will delete it!\n";
     	    close(OLIVIA);
     	    unlink "olivia.cutoff.$inputparameter.frequency.$frequency.$iframe"; 
     	    print "Delete Successful\n";
     	  }
     }


    }  ## if and else

    else {
     for $i (0..$natoms-1) {
      $line = <DATA>; chomp($line);
    	} 
    	}
    


     ####
  } ## if 
    $iframe++;
  } ## while

$totalframes=$iframe;
print "total number of frames is $totalframes\n";
print "total number of debris atom is $ndebris\n";
close(DATA);
}  ## foreach ends here.
