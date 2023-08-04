#!/usr/bin/perl

print"SUBROUTINE BUILD OK \n";

sub build{
	
	#local($moleculei)=@_;
	#use indice $moleculei and $nomtab
	
	$i2=$moleculei+1;
	if($nomtab eq "parent"){
		@getfgt=split('_',$mol[$moleculei]);
		@getfgt4=split(' ',$molping[$moleculei]);
	}
	else{
		@getfgt=split('_',$molchild[$moleculei]);
		@getfgt4=split(' ',$molpingchild[$moleculei]);
	};

		#BUILD	
		foreach $ani (0..@getfgt-1){
			print "molecule $nomtab $moleculei connection number $ani ($getfgt[$ani] from $mol[$moleculei] / $molping[$moleculei])\n" if($param{VERBOSITY} >= 1); 
			@getfgt2=split('-',$getfgt[$ani]); #$getfgt2[0] first fgt and $getfgt2[1] second one

			#first lego
				@getfgt3=split('\*',$getfgt2[0]); # in @getfgt3 0 = no fgt and in 1 = connecting atom
				$lelego3=$getfgt4[$getfgt3[0]-1]; # translate no fgt into no lego 
				$filelego=$lego[$lelego3];
				$filelegono=$legono[$lelego3];
				$anchor=$getfgt3[1];
				
				if($ani == 0){	
					$lego1="lego".$getfgt3[0];
					#print "$filelego $filelegono into $lego1\n";
					#chop($rep = `$leaexe/searchsdfbyno.pl $filelego $filelegono > $lego1` );
					&searchsdfif($filelego, $filelegono, $lego1);
					# clean file and take POINTS only 
					#chop($rep3 =`$leaexe/stripsdf $lego1 lego \"POINTS\"`);
					&stripsdf($lego1, "POINTS", "lego");
					rename "lego", $lego1;
				}
				else{ #next times take combin.sdf to continue the linking
				# get corresponding anchor
					$lego1="combin.sdf";
					open(LEG,"<$lego1");
					$flagname=0;
					$newanchor="";
					while(<LEG>){
						@get=split(' ',$_);
						if($flagname && $get[0] ne ''){
							#print "$_\n";
							if($get[1]==$anchor && $get[0]==$getfgt3[0]){
								$newanchor=$get[2];
							};	
						}
						else{
							$flagname=0;
						};	
						if ($_=~/^>/ && $_=~/<LEGO>/){
							$flagname=1;
						};
					};	
					#print "$anchor become $newanchor\n";
			       		$anchor=$newanchor;
				};
			
			#second lego	
			if($getfgt2[0] ne "" && $getfgt2[1] ne ""){
				@getfgt3=split('\*',$getfgt2[1]);
				$lelego2=$getfgt4[$getfgt3[0]-1]; # translate no fgt into no lego
				$filelego2=$lego[$lelego2];
				$filelegono2=$legono[$lelego2];
				$anchor2=$getfgt3[1];
				$lego2="lego".$getfgt3[0];
				#print "$filelego2 $filelegono2 into $lego2\n";
				#chop($rep2 = `$leaexe/searchsdfbyno.pl $filelego2 $filelegono2 > $lego2` );
				&searchsdfif($filelego2, $filelegono2, $lego2);
				# clean file and take POINTS only
				#chop($rep3 =`$leaexe/stripsdf $lego2 lego \"POINTS\"`);
				&stripsdf($lego2, "POINTS", "lego");
				rename "lego", $lego2;
				
				#link $lego1 and $lego2
				#print "link $lego1 $anchor $lego2 $anchor2\n";
				# DO NOT CHANGE THE ORDER lego 1 anchor lego2 anchor2
				chop($rep3 =`perl $leaexe/LINK_2MOL.pl $lego1 $anchor $lego2 $anchor2 2`);
				print "$rep3\n" if($param{VERBOSITY} >= 1);
				unlink $lego2;
				unlink $lego1 if($lego1 ne "combin.sdf");
			}
			elsif($getfgt2[0] ne "" && -e "$lego1"){ # moli is constituted by one fragment
				rename $lego1, "combin.sdf";
				last;
			};
		};
		#chop($rep3 =`$leaexe/stripsdf combin.sdf mol.sdf \"POINTS\"`); 
		&stripsdf("combin.sdf", "POINTS", "mol.sdf");
		unlink "combin.sdf";
		#chop($rep3 =`$leaexe/REPLACE_X.pl mol.sdf`);
		&replacex("mol.sdf","replacex.sdf");
		rename "replacex.sdf","mol.sdf";
		print "Create molecule no $i2 in mol.sdf DONE\n" if($param{VERBOSITY} >= 1);
	
};
