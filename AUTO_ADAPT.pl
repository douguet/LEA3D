#!/usr/bin/perl

print "SUBROUTINE AUTO_ADAPT.pl \n";

sub autoadapt{

###############################################################################

	#Adaptation of the selection pressure
	
		if($generation==0){
			$param{SCALEMIN}=2;
			#$param{SCALEMAX}=6;
			$param{SCALEMAX}=6;
			$scalemaxinit=$param{SCALEMAX} if($generation==0);
		};

		# Progressive pressure on best molecules: $param{SCALEMAX} is doubled at the end
		# weight for each molecule scale uses @scorepop
		
		$param{SCALEMAX} = $scalemaxinit*(1+(1-($param{GENMAX}-$generation)/$param{GENMAX}));

		print "\nSelection pressure [$param{SCALEMIN} ; $param{SCALEMAX}]\n";

###############################################################################

	#Adaptation of $param{PROBA_MUT} and $param{PROBA_CROSS} and $proba1+$proba2+$proba3+$proba4 
	#rq: if $sumfreq != 0 (in MAIN) then frequencies are taken into account by getlego subroutine (in MAIN) to randomly select a lego
	
		if($generation==0){
			$param{PROBA_MUT}=70;
			$param{PROBA_CROSS}=30; # In CHILD.pl, the probability is divided by 2 because it creates 2 children
			
			$proba1=25; #suppress
			$proba2=25; #add
			$proba3=25; #replace
			$proba4=25; #permutation

			# to allow child optimization
			$optimischild=0; # no child optimization for first part of the optimization
		}
		else{
		#$generation > 0
			# to allow child optimization
			if($generation > ($param{GENMAX}*2/3)){ #after two third 
				$optimischild=1;
			};
		};

		# permutation:
                if($optimischild){
                	$nbtestpermut=10;
                }
                else{
                	$nbtestpermut=1;
                };

		print "Mutation Types probability (percentage): \tsuppress one fgt/part= $proba1\tadd one fgt= $proba2\treplace one fgt= $proba3\tpermutation/move/scramble= $proba4\n";

		if($optimischild){
			print "Child optimization procedure selected\n";
		}
		else{
			print "Child optimization disabled, mutations are randomly applied\n";
		};


###############################################################################

	print "\n";
};
