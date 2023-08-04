#!/usr/bin/perl

#Keep the following print as the reply to the require perl function
print "";

sub cyclesdf{

########## FRAGMENTS (reconnaissance de cycles etc ....)

        $nbcycle=0;
        @cycle='';

	$debug=0;

# retrieve cycle

        foreach $cyc (1..$istratom){
		#print "depart $cyc\n";
                foreach $cyc3 (1..$istratom){
			$listpassage="";
			$listedirecte="";
			$contact1=-1;
			$suite="";
			$suitei="";
			$suitefinale="";
                        @gsuite4=split(' ',$ifonc[$cyc3]);
                        if(@gsuite4 > 1){
                                foreach $cyc2 (0..@gsuite4-1){
                                        $contact1=$cyc3 if($gsuite4[$cyc2] == $cyc);
                                };
                        };
			#$debug=1 if($cyc==6 && $contact1==5);
			#$debug=0 if($cyc==6 && $contact1!=5);
                	if($contact1 > -1){
				print "depart $cyc puis via $contact1\n" if($debug);
				@gsuite=split(' ',$ifonc[$contact1]);
                        	if(@gsuite > 1){
                                	foreach $cyc2 (0..@gsuite-1){
                                        	if($gsuite[$cyc2] != $cyc){
                                                	$listpassage=$listpassage." $gsuite[$cyc2] ";
                                                	$suite=$suite." $gsuite[$cyc2] ";
                                                	$suitei=$suitei." $contact1 ";
                                        	};
                                	};
                                	@gsuite2=split(' ',$suite);
                                	@gsuite3=split(' ',$suitei);
                                	$qs=@gsuite2-1;
                                	while($gsuite2[0] ne ""){
						print "\t suite $suite\n" if($debug);
						print "\t suitei $suitei\n" if($debug);
                                        	@fget=split(' ',$ifonc[$gsuite2[$qs]]);
						print "\t$ifonc[$gsuite2[$qs]]\n" if($debug);
                                        	$endsub=0;
                                        	foreach $cyc2 (0..@fget-1){
							if($listedirecte!~/ $gsuite2[$qs]-$fget[$cyc2] /){
							if($fget[$cyc2]==$cyc && $fget[$cyc2]!=$gsuite3[$qs] && $listpassage!~/ $fget[$cyc2] /){
								$endsub=1;
								print "$suite avant nulle\n" if($debug);
								print "$suitei avant nulle\n" if($debug);
								@gsuite5=split(' ',$suitei);
								foreach $qs3 (0..@gsuite5-1){
									if($suitefinale!~/ $gsuite5[@gsuite5-1-$qs3] /){
										$suitefinale=$suitefinale." $gsuite5[@gsuite5-1-$qs3] ";
									};
								};	
								$suitefinale=" $gsuite2[$qs] ".$suitefinale if($suitefinale!~/ $gsuite2[$qs] /);
								$suitefinale=$suitefinale." $cyc " if($suitefinale!~/ $cyc /);
								print "suite finale $suitefinale\n" if($debug);
								$cycle[$nbcycle]=$suitefinale;
								$nbcycle++;
								$suitefinale="";
								$listedirecte=$listedirecte." $gsuite2[$qs]-$fget[$cyc2] ";
							}	
                                                	elsif($fget[$cyc2]!=$gsuite3[$qs] && $listpassage!~/ $fget[$cyc2] / ){
                                                        	$endsub=1;
								print "\t add $fget[$cyc2]\n" if($debug);
                                                        	$listpassage=$listpassage." $fget[$cyc2] ";
                                                        	$suite=$suite." $fget[$cyc2] ";
                                                        	$suitei=$suitei." $gsuite2[$qs] ";
                                                	}
							else{
								print "\t no action for $fget[$cyc2] because = gsuite3[]\n" if($debug && $fget[$cyc2]==$gsuite3[$qs]);
								print "\t no action for $fget[$cyc2] because = listpassage ($listpassage)\n" if($debug && $listpassage=~/ $fget[$cyc2] /);
							};
							}
							else{
								print "search stop for $gsuite2[$qs]-$fget[$cyc2] because listedirecte=$listedirecte\n" if($debug);
							};
                                        	};
                                        	if($endsub==0){
                                                	$gsuite2[$qs]="";
                                                	$suite=join(' ',@gsuite2);
                                                	$suite=" ".$suite." ";
                                                	$gsuite3[$qs]="";
                                                	$suitei=join(' ',@gsuite3);
                                                	$suitei=" ".$suitei." ";
                                        	};
                                        	@gsuite2=split(' ',$suite);
                                        	@gsuite3=split(' ',$suitei);
                                        	$qs=@gsuite2-1;
                                	};
                        	};
                	};
		};
        };

#	foreach $lc (0..@cycle-1){
#	       print "$lc $cycle[$lc]\n";
#	};

###################################################################################
############### ELIMINATION DES DOUBLES

                $longc=$nbcycle-1;
                foreach $lc (0..$longc-2){
                        foreach $lc2 (($lc+1)..$longc-1){
                                if ($cycle[$lc2] ne ''){
                                        $nbzero=0;
                                        @get1= split(' ',$cycle[$lc]);
                                        @get2= split(' ',$cycle[$lc2]);
                                        if (@get1 == @get2){
                                                foreach $lc3 (0..@get1-1){
                                                        foreach $lc4 (0..@get2-1){
                                                                if ($get1[$lc3] == $get2[$lc4]){
                                                                        $get2[$lc4] = 0;
                                                                        $nbzero ++;
                                                                };
                                                        };
                                                };
                                                $cycle[$lc2]='' if ($nbzero == @get2);
                                        };
                                };
                        };
                };
#               print "\n";
                $nbc=0;
		$nbc56=0;
                @cyclef='';
                foreach $lc (0..$longc-1){
                        if ($cycle[$lc] ne ''){
				print "unique $cycle[$lc]\n" if($debug);
                                $cyclef[$nbc]=$cycle[$lc];
                                $nbc++;
				@get56=split(' ',$cycle[$lc]);
				$long56=@get56;
				$nbc56++ if($long56<=6);	
                        };
                };
		print "$nbc56 cycles a 3, 4, 5 ou 6 atomes\n" if($debug);
                print "$nbc cycles\n" if($debug);
                print "\n" if($debug);
};
