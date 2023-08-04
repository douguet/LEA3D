#!/usr/bin/perl

$leaexe=$0;
#windows requires 2 steps
$leaexe=~s/lea3d-CLASS_FGT\.pl$//;
$leaexe=~s/\/$//;
#print "perl scripts in $leaexe\n";

$key=$leaexe."/lea3d-KEY.pl";
require $key;

$fileadd='';
$f='';
$filekey='';
$optionx='';
$optionframe='';

my($f,$filekey,$optionx,$fileadd,$optionframe)=@ARGV;

if($f eq ''|| $filekey eq '' || $optionx eq '' || ($optionx ne "-" && $optionx ne "X" && $optionx ne "x") ){
	print "usage: classfgt <file.sdf> <file_key> <'X' or '-'> <output file.sdf> <option: '0' (default) or '1' for using the framework mode> \n";
	exit;
};

	$optionx="X" if($optionx eq "x");
     	&keymol2($f,$filekey,$optionx,$optionframe);
     	unlink "tmp" if(-e "tmp");
     	
	$fileadd="unique.sdf" if($fileadd eq '');
     	if(-e "diff.sdf" && !-z "diff.sdf"){
		#chop($nb = `$leaexe/NBSDF.pl diff.sdf`);
		$nb=&nbsdf("diff.sdf");
     		print "different molecules : $nb in $fileadd\n";
		rename "diff.sdf", "$fileadd";
	};	

	if(-e "same.sdf" && !-z "same.sdf"){
		#chop($nb = `$leaexe/NBSDF.pl same.sdf`);
		$nb=&nbsdf("same.sdf");
     		print "equivalentes molecules : $nb in same_$fileadd\n";
     		rename "same.sdf", "same_$fileadd";
	};	

     	if(-e "same_additional_key.sdf" && !-z "same_additional_key.sdf"){
		#chop($nb = `$leaexe/NBSDF.pl same_additional_key.sdf`);
		$nb=&nbsdf("same_additional_key.sdf");
     		print "equivalentes molecules to new generated keys : $nb in same_additional_key.sdf\n";
	};	
     
    	if(-e "exclude.sdf" && !-z "exclude.sdf"){ 
		#chop($nb = `$leaexe/NBSDF.pl exclude.sdf` );
		$nb=&nbsdf("exclude.sdf");
     		print "excluded molecules : $nb in exclude_$fileadd\n";
     		rename "exclude.sdf", "exclude_$fileadd";
	};	

###########################################################################################
#SUBROUTINES

sub nbsdf{
        local($fsdf)=@_;

        $nbligne=0;
        $nbd=0;
        open(IN,"<$fsdf");
        while(<IN>){
                $nbd++ if($_=~/\$\$\$\$/);
                $nbligne++;
        };
        close(IN);
        if($nbligne >=5 && $nbd==0){#case of .mol
                $nbd=1;
        };
        #print "$nbd\n";
        $nbd;
};

###################################################################################

