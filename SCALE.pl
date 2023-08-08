#!/usr/bin/perl

print"SUBROUTINE SCALE OK \n";

sub scale{

#	local($max,$min)=@_;

	$max=$param{SCALEMAX};
	$min=$param{SCALEMIN};

#***********************************************************************
#	UPDATE probabilities
#***********************************************************************

	$sumscore=0;
	$sumecart=0;
	$nbexactind=0;
	$fitmoymin=$scorepop[0];
	$fitmoymax=$scorepop[0];
	foreach $i (0..$param{POP}-1){
		$sumscore=$sumscore+$scorepop[$i];
		if($scorepop[$i] > $fitmoymax){
			$fitmoymax=$scorepop[$i];
		}
		elsif($scorepop[$i] < $fitmoymin){
			$fitmoymin=$scorepop[$i];
		};
		$nbexactind++;
	};
	$moyscore=$sumscore/($nbexactind);
	print"Mean/moyenne=$moyscore for $nbexactind\ molecules\n" if($param{VERBOSITY} >= 1);
	
	open(FIT,">>fitmoy.dat");
		#$sco1=sprintf"%3.2f",($moyscore/$sumw)*100;
		#$sco2=sprintf"%3.2f",($fitmoymax/$sumw)*100;
		#$sco3=sprintf"%3.2f",($fitmoymin/$sumw)*100;
		
		#@scorepop already in %
		$sco1=sprintf"%3.2f",$moyscore;
		$sco2=sprintf"%3.2f",$fitmoymax;
		$sco3=sprintf"%3.2f",$fitmoymin;
		printf FIT "$generation $sco1 $sco2 $sco3\n";
	close(FIT);

	foreach $i (0..$param{POP}-1){
		$moyb=($scorepop[$i]-$moyscore);
		$moyb=$moyb*$moyb;
		$sumecart=$sumecart+$moyb;
	};
	$ecarttype=sqrt($sumecart/$nbexactind);
	print"Standart error/Ecart-type=$ecarttype for $nbexactind molecules\n" if($param{VERBOSITY} >= 1);
	
	@poids="";
	$sumpoids=0;
	foreach $i (0..$param{POP}-1){
		$p=$scorepop[$i];
		if($ecarttype!=0){
			$p=($p-$moyscore)/$ecarttype;
		}
		else{
			$p=0;
		};
		$poidsp	=($p*(($max-$min)/2)+(($max+$min)/2));

		if($poidsp > 0){
			$poids[$i]=sprintf"%3.2f",$poidsp;
		}
		else{
			$poids[$i]=0.1;
			#print "neg $poidsp => poids $poids[$i]\n";
		};
		$sumpoids=$sumpoids+$poids[$i];
	};

	print "\n" if($param{VERBOSITY} >= 1);
	print "Weight:\n" if($param{VERBOSITY} >= 1);
	foreach $i (0..$param{POP}-1){
		print "mol[$i]\t$scorepop[$i]\=> $poids[$i]\n" if($param{VERBOSITY} >= 1);
	};
	print"\n" if($param{VERBOSITY} >= 1);
};
