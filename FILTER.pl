#!/usr/bin/perl


print"SUBROUTINE FILTER OK \n";

sub filter{

	local($file)=@_;	

	#chop($cms1 = `$leaexe/NBSDF.pl $file` );
	$cms1=&nbsdf($file);

	# apply FILTER
	if($param{FILTER}){
		print "FILTER by DRUG $paramdrug \n";
	
		chop($rep3 =`perl $leaexe/DRUG.pl $file $paramdrug 1`);
		#print "$rep3\n" if($param{VERBOSITY} >= 1);
		
		if(!-e "druglike.sdf"){
			print "Warning in filtering: all molecules are discarded !\n\n";
			unlink "notdruglike.sdf";
			unlink $file;
			$memo_remark[$memoi]=$memo_remark[$memoi]."filtering discards all molecules\n";
		}
		else{
			print "Remove non druglike\n\n" if($param{VERBOSITY} >= 1);
			rename ("druglike.sdf","$file");
			unlink "notdruglike.sdf";
			#chop($cms2 = `$leaexe/NBSDF.pl $file` );
			$cms2=&nbsdf($file);
			print "$file : $cms2 molecules after filtering\n" if($param{VERBOSITY} >= 1);
			if($cms2 != $cms1){
				print "$file : $cms1 molecules before filtering and $cms2 molecules after filtering\n";
				$memo_remark[$memoi]=$memo_remark[$memoi]."$cms1 molecules before filtering and $cms2 molecules after filtering\n";
			};	
		};			
	};
};

####################################################################################
