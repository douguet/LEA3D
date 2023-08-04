#!/usr/bin/perl

print"SUBROUTINE CHILD.pl OK \n";

sub child{

	# to allow child optimization: see AUTO_ADAPT.pl

	$sumoperator=$param{PROBA_CROSS}/2+$param{PROBA_MUT};

	#foreach $j ($param{ELITISM}..$param{POP}-1){
	$j=$param{ELITISM};
	while($j <= ($param{POP}-1)){
	# fill $molchild[]=$mol[] and $molpingchild[]=$molping[];

		$molchild[$j]="";
		$molpingchild[$j]="";
		
		$operator=int(rand($sumoperator))+1;
		
		$success=0;
		$nbessaimax=1000;
		
		#remark: $poids from mol0 mol1 ... not ordered !

# /////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////
 		
		#$param{PROBA_CROSS is divided by 2 to take account for both children
		if($param{POP} >= 2 && $operator <= ($param{PROBA_CROSS}/2) && $param{LONGGEN}!=1){
		#crossover
	
			#first parent
			$m1=-1;
			$nbessai=0; 
			while (($m1 == -1 || $mol[$m1]!~/-/) && $nbessai <= $nbessaimax){ # au moins une liaison !
				$nbessai++;
				$nb2=rand($sumpoids);
				$m1=0;
				$sumpar=$poids[0]; 
				while ($nb2>$sumpar){
					$m1++;
					$sumpar=$sumpar+$poids[$m1];
				};
			};	
			if($nbessai < $nbessaimax){
				#second parent
				$m2=-1;
				$nbessai2=0;
				while (($m2 == -1 || $mol[$m2]!~/-/ || $m2==$m1) && $nbessai2 <= $nbessaimax){
					$nbessai2++;
					$nb2=rand($sumpoids);
					$m2=0;
					$sumpar=$poids[0];
					while ($nb2>$sumpar){
						$m2++;
						$sumpar=$sumpar+$poids[$m2];
					};
				};
				if($nbessai2 < $nbessaimax && $mol[$m1]=~/-/ && $mol[$m2]=~/-/){
					print "\nparent1= $m1 and parent2= $m2    \n";
				
				# split parent 1	
					#dvt
					#$mol[$m1]="1*1-2*2_2*3-3*6_2*4-4*5_3*7-5*8_5*9-6*10";
					#$molping[$m1]="24 1 8 125 10 3";
					#$mol[$m1]="1*1-2*1_2*1-3*6_3*3-4*6_1*1-5*1";
					#$molping[$m1]=" 2 2 1 1  6 ";

					$fgt=$mol[$m1];
					$fgtping=$molping[$m1];
					@getfgt=split('_',$fgt);
					$break=int(rand(@getfgt));
					
					#dvt
					#$break=1;
					
					#print "\nmol$m1 $fgt / $fgtping\n";
					#print "break $break : $getfgt[$break]\n";
					
					&separate;		
					#@child1 new numbering scheme for 1st fragment
					#@child2 new numbering scheme for 2nd fragment
					#$childping1 and $childping2 no of legos
					
					@cross1="";
					@cross1b="";
					if($param{SCAFFOLD} ne '0' && $childping2=~/ 0 /){
						@cross1=@child2;
						@cross1b=@child1;
						$crossping1=$childping2;
						$crossping1b=$childping1;

						#$op_where[$j]="$getfgt4[$m1fgt2-1]";#from &separate
					}
					else{	
						@cross1=@child1;
						@cross1b=@child2; 
						$crossping1=$childping1;
						$crossping1b=$childping2;

						#$op_where[$j]="$getfgt4[$m1fgt1-1]";#from &separate
					};	

				# split parent 2	
				
					#dvt
					#$mol[$m2]="1*5-2*5_1*3-3*1_2*3-4*6_4*3-5*1";
					#$molping[$m2]=" 3 3  6  1 5";

					$fgt=$mol[$m2];
				       	$fgtping=$molping[$m2];	
					@getfgt=split('_',$fgt);
					$break=int(rand(@getfgt));
					#print "\nmol$m2 $fgt / $fgtping\n";
					#print "break $break : $getfgt[$break]\n";

					#dvt
					#$break=3;

					&separate;

					@cross2="";
					@cross2b="";
					if($param{SCAFFOLD} ne '0' && $childping2=~/ 0 /){
						@cross2=@child2;
						@cross2b=@child1;
						$crossping2=$childping2;
						$crossping2b=$childping1;

						#$op_where[$j]=$op_where[$j]."_$getfgt4[$m1fgt1-1]"; #from &separate
					}
					else{
						@cross2=@child1;
						@cross2b=@child2;
						$crossping2=$childping1;
						$crossping2b=$childping2;

						#$op_where[$j]=$op_where[$j]."_$getfgt4[$m1fgt2-1]"; #from &separate
					};	

				## First child
				# combine @cross1 and @cross2b (2 and 1b are not used so far)
					$longtab=@cross1;
				       	if($longtab==1){
						$molpingchild[$j]=$crossping1;
					}
					elsif($longtab > 1){	
						$molpingchild[$j]=$crossping1;
						foreach $k (1..@cross1-1){
							if($molchild[$j] eq ""){
								$molchild[$j]=$cross1[$k];
							}
							else{
								$molchild[$j]=$molchild[$j]."_".$cross1[$k];
							};	
						};	
					};
					#print "partial: molchild: $molchild[$j] / $molpingchild[$j]\n";
					@getfgt=split(' ',$molpingchild[$j]);
					$nbfgt=@getfgt;
				
					# connect 1 from cross1 to 1 from cross2b
					@getfgt2=split('\*',$cross2b[0]);
					$p1=$getfgt2[0]+$nbfgt;
					if($molchild[$j] eq ""){
						$molchild[$j]=$cross1[$k]."-".$p1."*$getfgt2[1]";
					}
					else{	
						$molchild[$j]=$molchild[$j]."_".$cross1[$k]."-".$p1."*$getfgt2[1]"; 
					};
					$molpingchild[$j]=$molpingchild[$j]." ".$crossping2b;
					
					foreach $k (1..@cross2b-1){
						@getfgt2=split('-',$cross2b[$k]);
						@getfgt3=split('\*',$getfgt2[0]);
						@getfgt5=split('\*',$getfgt2[1]);
						$p1=$getfgt3[0]+$nbfgt;
						$p2=$getfgt5[0]+$nbfgt;
						$molchild[$j]=$molchild[$j]."_".$p1."*$getfgt3[1]"."-".$p2."*$getfgt5[1]";
					};
					$molpingchild[$j]=~s/  / /g;
					$scorepopchild[$j]="not yet calculated";
					print "molchild[$j]: $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";

					$op_which[$j]="crossover";
					if($scorepop[$m1] >= $scorepop[$m2]){
						$op_scorep[$j]=$scorepop[$m1];
					}
					else{
					        $op_scorep[$j]=$scorepop[$m2];
					};

				## second child
				#combine @cross2 and @cross1b
					
					# next child
					$j++;
					$molchild[$j]="";
					$molpingchild[$j]="";

					$longtab=@cross2;
					if($longtab==1){
						$molpingchild[$j]=$crossping2;
					}
					elsif($longtab > 1){
						$molpingchild[$j]=$crossping2;
						foreach $k (1..@cross2-1){
							if($molchild[$j] eq ""){
								$molchild[$j]=$cross2[$k];
							}
							else{
								$molchild[$j]=$molchild[$j]."_".$cross2[$k];
							};
						};
					};
					#print "partial: molchild: $molchild[$j] / $molpingchild[$j]\n";
					
					@getfgt=split(' ',$molpingchild[$j]);
					$nbfgt=@getfgt;

					# connect 1 from cross2 to 1 from cross1b
					@getfgt2=split('\*',$cross1b[0]);
					$p1=$getfgt2[0]+$nbfgt;
					if($molchild[$j] eq ""){
						$molchild[$j]=$cross2[$k]."-".$p1."*$getfgt2[1]";
					}
					else{
						$molchild[$j]=$molchild[$j]."_".$cross2[$k]."-".$p1."*$getfgt2[1]";
					};	
					$molpingchild[$j]=$molpingchild[$j]." ".$crossping1b;

					foreach $k (1..@cross1b-1){
						@getfgt2=split('-',$cross1b[$k]);
						@getfgt3=split('\*',$getfgt2[0]);
						@getfgt5=split('\*',$getfgt2[1]);
						$p1=$getfgt3[0]+$nbfgt;
						$p2=$getfgt5[0]+$nbfgt;
						$molchild[$j]=$molchild[$j]."_".$p1."*$getfgt3[1]"."-".$p2."*$getfgt5[1]";
					};
					$molpingchild[$j]=~s/  / /g;
					$scorepopchild[$j]="not yet calculated";

					print "molchild[$j]: $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";

					$success=1;

					$op_which[$j]="crossover";
					if($scorepop[$m1] >= $scorepop[$m2]){
						$op_scorep[$j]=$scorepop[$m1];
					}
					else{
						$op_scorep[$j]=$scorepop[$m2];
					};
				}
				else{
					print "Select MUTATION: crossover failed in parent2 selection\n";
				};	
			}
			else{
				print "Select MUTATION: crossover failed in parent1 selection\n";
			};		
		};

# /////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////
 	
		$successm=0;
		if($operator > $param{PROBA_CROSS}/2 || $success==0){
			
			if($operator <= $param{PROBA_CROSS}/2){
				print "Select MUTATION: crossover not available Parameter POP = $param{POP} (must be >= 2)\n" if($param{POP} < 2);
				print "Select MUTATION: crossover not available Parameter LONGGEN = $param{LONGGEN} (must be 0 or > 1)\n" if($param{LONGGEN}==1);
			};
			
		#mutation
				
			#first parent	
			$nb2=rand($sumpoids);
			$m1=0;
			$sumpar=$poids[0];
			while ($nb2>$sumpar){
				$m1++;
				$sumpar=$sumpar+$poids[$m1];
			};
			#print "parent1= $m1\n";
			$fgt=$mol[$m1];
			$fgtping=$molping[$m1];
			@getsize=split(' ',$fgtping);
			$sizefgt=@getsize;

			#$probai can be modulated by autoadapt
			$typemutation=$proba1+$proba2+$proba3+$proba4;
			$quelmut=int(rand($typemutation))+1;
			if($quelmut <= $proba1){
				$nb2=1;
			}	
			elsif($quelmut <= ($proba1+$proba2)){
				$nb2=2;
			}
			elsif($quelmut <= ($proba1+$proba2+$proba3)){
				$nb2=3;
			}
			else{
				$nb2=4;
			};	
			if($nb2==2 && $param{LONGGEN}==$sizefgt && $param{LONGGEN} > 0){
				$nb2=3;
			};	
			
			#dvt
			#$nb2=2;
			
			#$nb2=1 suppress one fgt/part
			#$nb2=2 add one fgt
			#$nb2=3 replace one fgt
			#$nb2=4 either permutation of 2 substituents or move one substituent on another place
			#      or scramble around fgt with more than one <POINT>

####################################

			if($nb2==1){ #suppress one fragment, keep unchanged if nbfgt == 1
				print "\nsuppress fragment(s) into mol$m1 $fgt / $fgtping (score=$scorepop[$m1])\n";
				
				@getfgt=split('_',$fgt);
				$longtabsup=@getfgt;

				if($optimischild){

					#try all possibilities
					$bestscore=0;
					$bestmolpingchild="";
					$bestmolchild="";
					$besti=0;
					
					foreach $bri (0..$longtabsup-1){
						$break=$bri;

						&separate;

						@mut1="";
						@mut2="";
						$mutping1="";
						$mutping2="";
						@mut1=@child1;
						$mutping1=$childping1;
						@mut2=@child2;
						$mutping2=$childping2;

						foreach $brj (1..2){
							$molchild[$j]="";
							$molpingchild[$j]="";
						
							if($brj==1){
								@mut="";
								@mut=@mut1;
								$mutping=$childping1;

								#$op_where[$j]="$getfgt4[$m1fgt1-1]"; # from &separate
							}
							else{
								@mut="";
								@mut=@mut2;
								$mutping=$childping2;
								
								#$op_where[$j]="$getfgt4[$m1fgt2-1]"; # from &separate
							};

							if($param{SCAFFOLD} eq '0' || ($param{SCAFFOLD} ne '0' && $mutping=~/ 0 /)){
								$longtab=@mut;
								if($longtab==1){
									$molpingchild[$j]=$mutping;
									$molchild[$j]=$mut[0];
								}
								else{
									$molpingchild[$j]=$mutping;
									foreach $k (1..@mut-1){
										if($molchild[$j] eq ""){
											$molchild[$j]=$mut[$k];
										}		
										else{
											$molchild[$j]=$molchild[$j]."_".$mut[$k];
										};
									};	
								};		
								$molpingchild[$j]=~s/  / /g;

								#Build, Optimis and Score this child
								print "test $molchild[$j] / $molpingchild[$j]\n";
								&bos($j,"child");
								if(-e "mol.sdf" && !-z "mol.sdf"){
									$besti++;
									print "child besti $besti score=$score[$ranking[1]]\n";
									if($bestscore < $score[$ranking[1]]){
										$bestscore=$score[$ranking[1]];
										$bestmolchild=$molchild[$j];
										$bestmolpingchild=$molpingchild[$j];
									};
								};	
							};	
						};#foreach $brj	

					};#foreach $bri
				};

				if($optimischild && $bestmolchild ne ""){
					$successm=1;
					$molchild[$j]=$bestmolchild;
					$molpingchild[$j]=$bestmolpingchild;
					$scorepopchild[$j]=$bestscore;
					print "Child $j ($besti tests): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
				}
				elsif($optimischild==0 || $bestmolchild eq ""){
					# random
					
					$molchild[$j]="";	
					$molpingchild[$j]="";

					@getfgt=split('_',$fgt); 
					$break=int(rand(@getfgt));
					&separate;
					print "keep the more complex => child1: @child1 / $childping1 ; child2: @child2 / $childping2\n";
					$longtab=@child1;
					$longtab2=@child2;
					if(($param{SCAFFOLD} eq '0' && $longtab >= $longtab2) || ($param{SCAFFOLD} ne '0' && $childping1=~/ 0 /)){
						@mut="";
						@mut=@child1; #keep first part
						$mutping=$childping1;

						 #$op_where[$j]="$getfgt4[$m1fgt1-1]"; # from &separate
					}
					else{
						@mut="";
						@mut=@child2; #keep second part
						$mutping=$childping2;

						#$op_where[$j]="$getfgt4[$m1fgt2-1]"; # from &separate
					};

					$longtab=@mut;
					if($longtab==1){
						$molpingchild[$j]=$mutping;
						$molchild[$j]=$mut[0];
					}
					else{
						$molpingchild[$j]=$mutping;
						foreach $k (1..@mut-1){
							if($molchild[$j] eq ""){
								$molchild[$j]=$mut[$k];
							}
							else{
								$molchild[$j]=$molchild[$j]."_".$mut[$k];
							};
						};

					};	
					$molpingchild[$j]=~s/  / /g;
					$scorepopchild[$j]="not yet calculated";
							
					if($molchild[$j] eq ""){
						print "\nWarning! suppress one fragment in mol$m1 failed !\n";
					}
					else{
						$successm=1;
						print "molchild[$j] (random): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
					};
				};	
			
				#STUDIES on operator efficiency
				$op_which[$j]="suppress";
				$op_scorep[$j]=$scorepop[$m1];
				if($successm){
					$op_where[$j]="";
					@brique=split(' ',$fgtping);
					$testbrique=" ".$molpingchild[$j]." ";
					foreach $briquei (0..@brique-1){
						# which lego has been removed
						if($testbrique!~/ $brique[$briquei] /){
							$op_where[$j]=$op_where[$j]." ".$brique[$briquei];
						};
					};
				};

####################################
			}
			elsif($nb2==2){ #add one fragment as in initialize pop in MAIN. if no free anchor: keep unchanged

				$molchild[$j]="";
				$molpingchild[$j]="";

				#try all possibilities
				$bestscore=0;
				$bestmolpingchild="";
				$bestmolchild="";
				$besti=0;

				#preliminary test
				$tmpa=&getfreeanchor($fgt,$fgtping); #null if no free anchor !

				if($tmpa eq ""){
					print "\nadd one fragment into mol$m1 $fgt / $fgtping failed: no free anchor\n";
				}	
				else{	
					print "\nadd one fragment into mol$m1 $fgt / $fgtping (score=$scorepop[$m1])\n";
					@getfgt=split(' ',$fgtping);
					$newfgt=@getfgt;
					$newfgt++;

					if($optimischild){
						#dvt
						#print "Add: subset created from privileged legos\n";
						
						$tmpmol=&getlego("any");
						$tmpmol=~s/(.*)\*(.*)/$1/;
						$tmpb=$newfgt."*".$2;
						print "try and score lego no $tmpmol at anchor position $2 -> fgt number $newfgt in molchild$j\n";
						
						#try all free anchor
						$listtmpa=&getfreeanchor($fgt,$fgtping,"list");
						@getlistfreeanchor=split(' ',$listtmpa);
						foreach $brj (0..@getlistfreeanchor-1){

							#new lego each time ?
							#$tmpmol=&getlego("any");
							#$tmpmol=~s/(.*)\*(.*)/$1/;
							#$tmpb=$newfgt."*".$2;
							#print "try and score lego no $tmpmol at anchor position $2 -> fgt number $newfgt in molchild$j\n";
							
							if($newfgt==2){
								$molchild[$j]=$getlistfreeanchor[$brj]."-".$tmpb;
								$molpingchild[$j]=$fgtping." $tmpmol";
							}
							else{
								$molchild[$j]=$fgt."_".$getlistfreeanchor[$brj]."-".$tmpb;
								$molpingchild[$j]=$fgtping." $tmpmol";
							};
							#Build, Optimis and Score this child
							print "$molchild[$j] / $molpingchild[$j]\n";
							&bos($j,"child");
							if(-e "mol.sdf" && !-z "mol.sdf"){
								$besti++;
								print "child besti $besti score=$score[$ranking[1]]\n";
								if($bestscore < $score[$ranking[1]]){
									$bestscore=$score[$ranking[1]];
									$bestmolchild=$molchild[$j];
									$bestmolpingchild=$molpingchild[$j];
								};
							};
						};
					};#$optimischild
					
					if($optimischild && $bestmolchild ne ""){
						$successm=1;
						$molchild[$j]=$bestmolchild;
						$molpingchild[$j]=$bestmolpingchild;
						$scorepopchild[$j]=$bestscore;
						print "Child $j ($besti tests): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
					}
					elsif($optimischild==0 || $bestmolchild eq ""){
						# random
						
						$molchild[$j]="";
						$molpingchild[$j]="";

						#print "no lego bring >= $pcvumax \%: try random selection\n" if($bestmolchild eq "" && $optimischild==1 && $nbdiff!=0);
						#print "difference in GC = 0 -> try random selection\n" if($nbdiff==0 && $optimischild==1);

						$tmpmol=&getlego("any");
						$tmpmol=~s/(.*)\*(.*)/$1/;
						$tmpb=$newfgt."*".$2;
						print "select lego $tmpmol at anchor position $2 (from anchors $points[$tmpmol])-> fgt number $newfgt in molchild[$j]\n";
						 
						if($newfgt==2){
							$molchild[$j]=$fgt."-".$tmpb;
							$molpingchild[$j]=$fgtping." $tmpmol";
						}
						else{
							$molchild[$j]=$fgt."_".$tmpa."-".$tmpb;
							$molpingchild[$j]=$fgtping." $tmpmol";
						};
						$scorepopchild[$j]="not yet calculated";

						if($molchild[$j] eq ""){
							print "\nWarning! add one fragment into mol$m1 $fgt / $fgtping failed\n";
						}
						else{
							$successm=1;
							print "molchild[$j] (random): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
						};

					};

				}; #if($tmpa eq "")

				#STUDIES on operator efficiency
				$op_which[$j]="add";
				$op_scorep[$j]=$scorepop[$m1];
				if($successm){
					$op_where[$j]="";
					@brique=split(' ',$molpingchild[$j]);
					$testbrique=" ".$fgtping." ";
					foreach $briquei (0..@brique-1){
						# which lego has been added:
						if($testbrique!~/ $brique[$briquei] /){
							$op_where[$j]=$op_where[$j]." ".$brique[$briquei];
						};
					};
				};	


####################################				
			}
			elsif($nb2==3){ #replace one fragment
				print "\nreplace one fragment into mol$m1 $fgt / $fgtping (score=$scorepop[$m1])\n";
			
				$molchild[$j]="";
				$molpingchild[$j]="";

				#try all possibilities
				$bestscore=0;
				$bestmolchild="";
				$bestmolpingchild="";
				$besti=0;
			
				if($optimischild){
					
					@getfgt=split(' ',$fgtping);
					$longtabreplace=@getfgt;
					$brmlower=0;
					foreach $brm ($brmlower..$longtabreplace-1){

						$molchild[$j]="";
						$molpingchild[$j]="";

						$selectfgt=$brm; #position $selectfgt and fgt no ($selectfgt+1)
						@getfgt=split(' ',$fgtping);
						$selectnolego=$getfgt[$selectfgt];#lego no $selectnolego
						$selectfgtinmol=$selectfgt+1; #fragment number in mol[]

						if($param{SCAFFOLD} ne '0' && $selectnolego==0){
							print "replace this fragment into mol$m1 $fgt / $fgtping failed, scaffold lego 0 cannot be replaced !\n";
						}	
						else{

							#$op_where[$j]="$selectnolego";

							#constraint = number of occupied anchors
							$whichfgtbis=" ".$fgt." ";
							$whichfgtbis=~s/-/ /g;
							$whichfgtbis=~s/_/ /g;
							$car=$whichfgtbis;
							$selectfgt2=$selectfgt+1;
							$nocc=$car=~s/ $selectfgt2\*/ $selectfgt2\*/g;
							$constraint=$nocc; # constraint on anchors
					
							print "replace fgt no $selectfgtinmol (lego $selectnolego, $typelego[$selectnolego] anchors and $constraint substituted) with new fgt of the minimal number of $constraint anchors\n";

							if($param{LEGODB}){
                                                        	print "get similar lego $selectfgtinmol, $fgt, $fgtping\n";
                                                                $tmpmol=&getsimilarlego($selectfgtinmol,$fgt,$fgtping);
                                                        }
                                                        else{
                                                                $tmpmol=&getlego($constraint);
                                                                $tmpmol=~s/(.*)\*(.*)/$1/;
                                                        };
                                                	if($tmpmol ne ""){
                                                        	print "select new lego $tmpmol ($typelego[$tmpmol] anchors: $points[$tmpmol])\n";
                                                        	@getfgt=split(' ',$fgtping);
                                                        	$getfgt[$selectfgt]=$tmpmol;
                                                        	$molpingchild[$j]=join(' ',@getfgt);

                                                        	@listpoint=split(' ',$points[$tmpmol]);
                                                        	@getfgt=split('_',$fgt);
                                                        	foreach $k (0..@getfgt-1){
                                                                	@getfgt2=split('-',$getfgt[$k]);
                                                                	foreach $l (0..@getfgt2-1){
                                                                        	@getfgt3=split('\*',$getfgt2[$l]);
                                                                        	if($getfgt3[0] == $selectfgtinmol){
                                                                                	$vu=-1;
                                                                                	foreach $m (0..@listpoint-1){
                                                                                        	if($vu==-1 && $listpoint[$m] ne ""){
                                                                                                	$vu=$listpoint[$m];
                                                                                                	$listpoint[$m]="";
                                                                                        	};
                                                                                	};
                                                                                	$getfgt3[1]=$vu;
                                                                                	print "warning inconsistencies with anchors ?\n" if($vu==-1);
                                                                        	};
                                                                        	$getfgt2[$l]=$getfgt3[0]."*".$getfgt3[1];
                                                                	};
                                                                	$getfgt[$k]=join('-',@getfgt2);
                                                        	};
                                                        	$molchild[$j]=join('_',@getfgt);

								#Build, Optimis and Score this child
								print "$molchild[$j] / $molpingchild[$j]\n";
								&bos($j,"child");
								if(-e "mol.sdf" && !-z "mol.sdf"){
									$besti++;
									print "child besti $besti score=$score[$ranking[1]]\n";
									if($bestscore < $score[$ranking[1]]){
										$bestscore=$score[$ranking[1]];
										$bestmolchild=$molchild[$j];
										$bestmolpingchild=$molpingchild[$j];
									};
								};

                                                	};#if($tmpmol ne "")
						};# if($param{SCAFFOLD} ne '0' && $selectnolego==0)
					}; #foreach $brm	
				};#if($optimischild)

				if($optimischild && $bestmolchild ne ""){
					$successm=1;
					$molchild[$j]=$bestmolchild;
					$molpingchild[$j]=$bestmolpingchild;
					$scorepopchild[$j]=$bestscore;
					print "Child $j ($besti tests): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
				}
				elsif($optimischild==0 || $bestmolchild eq ""){
					# random
				
					$molchild[$j]="";
					$molpingchild[$j]="";

					#print "no lego bring >= $pcvumax \%: try random selection\n" if($bestmolchild eq "" && $optimischild==1 && $nbdiff!=0);
					#print "difference in GC = 0 -> try random selection\n" if($nbdiff==0 && $optimischild==1);

					@getfgt=split(' ',$fgtping);
					$longtab=@getfgt;
					$selectfgt=int(rand(@getfgt)); #position $selectfgt and fgt no ($selectfgt+1)
					$selectnolego=$getfgt[$selectfgt];
					$selectfgtinmol=$selectfgt+1; #fragment number in mol[]
					
					if($param{SCAFFOLD} ne '0'){
						if($longtab > 1){
							while($selectnolego==0){
								$selectfgt=int(rand(@getfgt));
								$selectnolego=$getfgt[$selectfgt];
								$selectfgtinmol=$selectfgt+1; #fragment number in mol[]
							};
						}
						else{
							print "replace one fragment into mol$m1 $fgt / $fgtping failed, scaffold lego 0 cannot be replaced !\n";
						};
					};
					if($selectnolego != 0){
						$whichfgtbis=" ".$fgt." ";
						$whichfgtbis=~s/-/ /g;
						$whichfgtbis=~s/_/ /g;
						$car=$whichfgtbis;
						$selectfgt2=$selectfgt+1;
						$nocc=$car=~s/ $selectfgt2\*/ $selectfgt2\*/g;
						$constraint=$nocc;

						print "replace fgt no $selectfgtinmol (lego $selectnolego, $typelego[$selectnolego] anchors and $constraint substituted) with new fgt of the minimal number of $constraint anchors\n";

						$tmpmol=&getlego($constraint);
						$tmpmol=~s/(.*)\*(.*)/$1/;
						if($tmpmol ne ""){
							print "select new lego $tmpmol ($typelego[$tmpmol] anchors: $points[$tmpmol])\n";
							@getfgt=split(' ',$fgtping);
							$getfgt[$selectfgt]=$tmpmol;
							$molpingchild[$j]=join(' ',@getfgt);
				
							@listpoint=split(' ',$points[$tmpmol]);
							@getfgt=split('_',$fgt);	
							foreach $k (0..@getfgt-1){
								@getfgt2=split('-',$getfgt[$k]);
								foreach $l (0..@getfgt2-1){
									@getfgt3=split('\*',$getfgt2[$l]);
									if($getfgt3[0] == $selectfgtinmol){
										$vu=-1;
										foreach $m (0..@listpoint-1){
											if($vu==-1 && $listpoint[$m] ne ""){
												$vu=$listpoint[$m];
												$listpoint[$m]="";
											};	
										};
										$getfgt3[1]=$vu;
										print "warning inconsistencies with anchors ?\n" if($vu==-1);
									};
									$getfgt2[$l]=$getfgt3[0]."*".$getfgt3[1];
								};
								$getfgt[$k]=join('-',@getfgt2);
							};
							$molchild[$j]=join('_',@getfgt);
							$scorepopchild[$j]="not yet calculated";

						};#if($tmpmol ne "")

					};#if($selectnolego != 0){

					if($molchild[$j] eq ""){
							print "\nWarning! replace one fragment into mol$m1 $fgt / $fgtping failed\n";
					}
					else{
						$successm=1;
						print "molchild[$j] (random): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
					};	

				};

				#STUDIES on operator efficiency
				$op_which[$j]="replace";
				$op_scorep[$j]=$scorepop[$m1];
				if($successm){
					$op_where[$j]="";
					@brique=split(' ',$molpingchild[$j]);
					$testbrique=" ".$fgtping." ";
					foreach $briquei (0..@brique-1){
					# which lego has been added:
						if($testbrique!~/ $brique[$briquei] /){
							$op_where[$j]=$op_where[$j]." ".$brique[$briquei];
						};
					};	
				};


####################################				
			}	
			elsif($nb2==4){ #permutation of one or several fragment
				print "\npermutation of fragments into mol$m1 $fgt / $fgtping (score=$scorepop[$m1])\n";

				#try $nbtestpermut possibilities
				$bestscore=0;
				$bestmolpingchild="";
				$bestmolchild="";
				$besti=0;

				$molchild[$j]="";
				$molpingchild[$j]="";

				@getfgt=split(' ',$fgtping);
				$longtab=@getfgt;

				if($longtab > 1){
				foreach $permuti (1..$nbtestpermut){

					@getfgt=split(' ',$fgtping);
					$selectfgt=int(rand(@getfgt)); #position $selectfgt and fgt no ($selectfgt+1)
					$selectnolego=$getfgt[$selectfgt];
					$constraint=$typelego[$selectnolego];
					$selectfgt++; # no of fgt in mol$i

					print "change the position of fgt no $selectfgt = lego $selectnolego (anchors= $points[$selectnolego])\n";
				
					if($constraint==1 && $longtab > 1){
				        	#print "find another fgt with only one substituant\n";	
					
						# try to find another fgt with only one substituant
						$whichfgtbis=" ".$fgt." ";
				       		$whichfgtbis=~s/-/ /g;
						$whichfgtbis=~s/_/ /g;
					
						$forceend=0;
						$vu=-1;
				       		while($forceend < 10000 && $vu==-1){
							$selectfgt2=int(rand(@getfgt));
							$selectfgt2++;
							$car=$whichfgtbis;
							$nocc=$car=~s/ $selectfgt2\*/ $selectfgt2\*/g;
							if($nocc==1 && $selectfgt2 != $selectfgt){
								$vu=$selectfgt2-1; # $vu = indice
							};
							$forceend++;	
						};
						#dvt
						#$vu=-1;
						if($vu != -1){
							@listpoint=split(' ',$points[$selectnolego]); #constraint==1 then $listpoint[0] is the anchor no
							$selectnolego2=$getfgt[$vu];
							$selectfgt2=$vu+1;	
							@getfgt=split('_',$fgt);
							foreach $k (0..@getfgt-1){
								@getfgt2=split('-',$getfgt[$k]);
								foreach $l (0..@getfgt2-1){
									@getfgt3=split('\*',$getfgt2[$l]);
									$anchor2=$getfgt3[1] if($getfgt3[0]== $selectfgt2);
								};
							};	
							print "permutation with terminal fgt no $selectfgt2 anchor $anchor2\n";	
							foreach $k (0..@getfgt-1){
								@getfgt2=split('-',$getfgt[$k]);
								foreach $l (0..@getfgt2-1){
									@getfgt3=split('\*',$getfgt2[$l]);
									if($getfgt3[0]== $selectfgt){
										$getfgt3[0]=$selectfgt2;
										$getfgt3[1]=$anchor2;
									}	
									elsif($getfgt3[0]==$selectfgt2){
										$getfgt3[0]=$selectfgt;
										$getfgt3[1]=$listpoint[0];
									};
									$getfgt2[$l]=$getfgt3[0]."*".$getfgt3[1];
								};
								$getfgt[$k]=join('-',@getfgt2);
							};	
							$molchild[$j]=join('_',@getfgt);
							$molpingchild[$j]=$fgtping;

							if($optimischild){
								#Build, Optimis and Score this child
								print "test $molchild[$j] / $molpingchild[$j]\n";
								&bos($j,"child");
								if(-e "mol.sdf" && !-z "mol.sdf"){
									$besti++;
									print "child besti $besti score=$score[$ranking[1]]\n";
									if($bestscore < $score[$ranking[1]]){
										$bestscore=$score[$ranking[1]];
										$bestmolchild=$molchild[$j];
										$bestmolpingchild=$molpingchild[$j];
									};
								};
							}
							else{
								$scorepopchild[$j]="not yet calculated";
							};
							#$op_where[$j]=$selectnolego."_".$selectnolego2;
						}
						else{
							#failed to find another fgt with only one substituant thus search for a new anchor in mol$i where to bond $selectnolego
						
							print "Cannot find fgt with only one substituant: search a new anchor point\n";
							$tmpa=&getfreeanchor($fgt,$fgtping);
							$tmpa=~s/(.*)\*(.*)/$1/;
							$tmpanchor=$2;
							if($tmpa eq "" || $tmpa == $selectfgt){ #null or only 1 fgt
								print "no free anchor, moving of fgt failed ($tmpa)\n";
							}
							else{	
								print "free anchor in mol$m1 $tmpa (anchor $tmpanchor)\n";
								@listpoint=split(' ',$points[$selectnolego]);#$listpoint[0] is the anchor no
								@getfgt=split('_',$fgt);
								foreach $k (0..@getfgt-1){
									$supp=0;
									@getfgt2=split('-',$getfgt[$k]);
									foreach $l (0..@getfgt2-1){
										@getfgt3=split('\*',$getfgt2[$l]);
										if($getfgt3[0]== $selectfgt){
											$supp=1;
										};	
									};
									$getfgt[$k]="" if($supp);
								};
								$molchild[$j]=join('_',@getfgt);
								$molchild[$j]=~s/^_//;
								if($molchild[$j] eq ""){
									$molchild[$j]=$tmpa."*".$tmpanchor."-".$selectfgt."*".$listpoint[0];
								}
								else{	
									$molchild[$j]=$molchild[$j]."_".$tmpa."*".$tmpanchor."-".$selectfgt."*".$listpoint[0];
								};
								$molpingchild[$j]=$fgtping;	

								if($optimischild){
									#Build, Optimis and Score this child
									print "test $molchild[$j] / $molpingchild[$j]\n";
									&bos($j,"child");
									if(-e "mol.sdf" && !-z "mol.sdf"){
										$besti++;
										print "child besti $besti score=$score[$ranking[1]]\n";
										if($bestscore < $score[$ranking[1]]){
											$bestscore=$score[$ranking[1]];
											$bestmolchild=$molchild[$j];
											$bestmolpingchild=$molpingchild[$j];
										};
									};
								}
                                                        	else{
                                                                	$scorepopchild[$j]="not yet calculated";
								};
								#$op_where[$j]=$selectnolego;
							};
						};	
					}#if($constraint==1 && $longtab > 1){
					else{# scramble substitution pattern
						print "scramble around fragment $selectfgt\n";
						@listpoint=split(' ',$points[$selectnolego]);
						@getfgt=split('_',$fgt);
						foreach $k (0..@getfgt-1){
							@getfgt2=split('-',$getfgt[$k]);
							foreach $l (0..@getfgt2-1){
								@getfgt3=split('\*',$getfgt2[$l]);
								if($getfgt3[0] == $selectfgt){
									$vu=-1;
									while($vu==-1){
										$nb3=int(rand($constraint));
										if($listpoint[$nb3] ne ""){
											$vu=$listpoint[$nb3];
											$listpoint[$nb3]="";
										};	
									};
									$getfgt3[1]=$vu;
									#print "select randomly anchor $vu\n";
									print "warning inconsistencies with anchors ?\n" if($vu==-1);	
								};
								$getfgt2[$l]=$getfgt3[0]."*".$getfgt3[1];
							};
							$getfgt[$k]=join('-',@getfgt2);
						};
						$molchild[$j]=join('_',@getfgt);
						$molpingchild[$j]=$fgtping;

						if($optimischild){
							#Build, Optimis and Score this child
							print "test $molchild[$j] / $molpingchild[$j]\n";
							&bos($j,"child");
							if(-e "mol.sdf" && !-z "mol.sdf"){
								$besti++;
								print "child besti $besti score=$score[$ranking[1]]\n";
								if($bestscore < $score[$ranking[1]]){
									$bestscore=$score[$ranking[1]];
									$bestmolchild=$molchild[$j];
									$bestmolpingchild=$molpingchild[$j];
								};
							};
						}
                                                else{
                                                	$scorepopchild[$j]="not yet calculated";	
						};
						#$op_where[$j]=$selectnolego;
					};

					$op_which[$j]="permutation";
					$op_scorep[$j]=$scorepop[$m1];
					
				};#foreach $permuti
				};#if longtab > 1

				if($optimischild && $bestmolchild ne ""){
					$successm=1;
					$molchild[$j]=$bestmolchild;
					$molpingchild[$j]=$bestmolpingchild;
					$scorepopchild[$j]=$bestscore;
					print "molchild[$j] Child $j ($besti tests): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n" if($optimischild);
				}
				elsif($optimischild==0 && $molchild[$j] ne ""){
					$successm=1;
					print "molchild[$j] Child $j (random): select $molchild[$j] / $molpingchild[$j] (score=$scorepopchild[$j])\n";
				}
				elsif($optimischild && $bestmolchild eq ""){
					print "\nWarning! permutation in mol$m1 failed !\n";
				};


			};#end of $nb=4	

		};#end of proba_mutation

		@getsize=split(' ',$molpingchild[$j]);
		$sizefgt=@getsize;
		
		if($param{LONGGEN} > 0){
			#$j++ if(($success || $successm) && $sizefgt > $param{LONGGEN});
			$j++ if(($success || $successm) && $sizefgt <= $param{LONGGEN});#Error corrected ? Feb2017
		}
		else{
			$j++ if($success || $successm);
		};	
	};
	
};

##############################################################################################################

sub separate {

	print "$fgt / $fgtping (break at $break)\n";
	
	@getfgt=split('_',$fgt);
	@getfgt4=split(' ',$fgtping);
	
	@getfgt2=split('-',$getfgt[$break]); #get both fragments
	@getfgt3=split('\*',$getfgt2[0]);#get 1st fragment number $getfgt3[0]
	@getfgt5=split('\*',$getfgt2[1]);#get 2nd fragment number $getfgt5[0]	
			
	$m1fgt1=$getfgt3[0];
	$m1fgt2=$getfgt5[0];
					
	@part1=""; #keep fragment number
	@part2=""; #discard fragment number
	$part1[0]=$m1fgt1;
	$part2[0]=$m1fgt2;
					
	@child1="";
	$child1[0]="1*".$getfgt3[1];
	$childping1=" ".$getfgt4[$m1fgt1-1];
	@child2="";
	$child2[0]="1*".$getfgt5[1];
	$childping2=" ".$getfgt4[$m1fgt2-1];
		
	#split bounds
	$part1i=1;
	$part2i=1;
	$nbfgt=@getfgt4;

	@vugetfgt="";
	$allvu=0;
	$countvu=0;
	while($allvu < $nbfgt && $countvu <100){
	foreach $j2 (0..@getfgt-1){
		if($j2 != $break && $vugetfgt[$j2] eq ""){
			@getfgt2=split('-',$getfgt[$j2]);
			@getfgt3=split('\*',$getfgt2[0]);
			@getfgt5=split('\*',$getfgt2[1]);
			if($getfgt3[0]==$m1fgt1){
				$part1[$part1i]=$getfgt5[0];
				$p=$part1i+1;
				$child1[$part1i]="1*"."$getfgt3[1]"."-"."$p"."*"."$getfgt5[1]";
				$childping1=$childping1." $getfgt4[$getfgt5[0]-1]";
				$part1i++;
				$vugetfgt[$j2]=1;
			}
			elsif($getfgt5[0]==$m1fgt1){
				$part1[$part1i]=$getfgt3[0];
				$p=$part1i+1;
				$child1[$part1i]="1*"."$getfgt5[1]"."-"."$p"."*"."$getfgt3[1]";
				$childping1=$childping1." $getfgt4[$getfgt3[0]-1]";
				$part1i++;
				$vugetfgt[$j2]=1;
			}
			elsif($getfgt3[0]==$m1fgt2){
				$part2[$part2i]=$getfgt5[0];
				$p=$part2i+1;
				$child2[$part2i]="1*"."$getfgt3[1]"."-"."$p"."*"."$getfgt5[1]";
				$childping2=$childping2." $getfgt4[$getfgt5[0]-1]";
				$part2i++;
				$vugetfgt[$j2]=1;
			}
			elsif($getfgt5[0]==$m1fgt2){
				$part2[$part2i]=$getfgt3[0];
				$p=$part2i+1;
				$child2[$part2i]="1*"."$getfgt5[1]"."-"."$p"."*"."$getfgt3[1]";
				$childping2=$childping2." $getfgt4[$getfgt3[0]-1]";
				$part2i++;
				$vugetfgt[$j2]=1;
			}
			else{ #neither $m1fgt1 nor $m1fgt2 but necessary another $part1 or $part2
				$vu=0;
				foreach $k (0..@part1-1){
					if($part1[$k]==$getfgt3[0]){
						$part1[$part1i]=$getfgt5[0];
						$p=$part1i+1;
						$k2=$k+1;
						$child1[$part1i]="$k2"."*"."$getfgt3[1]"."-"."$p"."*"."$getfgt5[1]";
						$childping1=$childping1." $getfgt4[$getfgt5[0]-1]";
						$part1i++;
						$vugetfgt[$j2]=1;
						$vu=1;
					}
					elsif($part1[$k]==$getfgt5[0]){
						$part1[$part1i]=$getfgt3[0];
						$p=$part1i+1;
						$k2=$k+1;
						$child1[$part1i]="$k2"."*"."$getfgt5[1]"."-"."$p"."*"."$getfgt3[1]";
						$childping1=$childping1." $getfgt4[$getfgt3[0]-1]";
						$part1i++;
						$vugetfgt[$j2]=1;
						$vu=1;
					};
				};	
				if($vu==0){
					foreach $k (0..@part2-1){
						if($part2[$k]==$getfgt3[0]){
							$part2[$part2i]=$getfgt5[0];
							$p=$part2i+1;
							$k2=$k+1;
							$child2[$part2i]="$k2"."*"."$getfgt3[1]"."-"."$p"."*"."$getfgt5[1]";
							$childping2=$childping2." $getfgt4[$getfgt5[0]-1]";
							$part2i++;
							$vugetfgt[$j2]=1;
							$vu=1;
						}
						elsif($part2[$k]==$getfgt5[0]){
							$part2[$part2i]=$getfgt3[0];
							$p=$part2i+1;
							$k2=$k+1;
							$child2[$part2i]="$k2"."*"."$getfgt5[1]"."-"."$p"."*"."$getfgt3[1]";
							$childping2=$childping2." $getfgt4[$getfgt3[0]-1]";
							$part2i++;
							$vu=1;	
							$vugetfgt[$j2]=1;
						};
					};
				};
			};	
			$allvu=@part2+@part1;
		};		
		$countvu++;
	};
	};
	
	$childping1=$childping1." ";
	$childping2=$childping2." ";

	$allvu=@part2+@part1;
	print "warning did not correctly separate fgts\n" if($allvu < $nbfgt);
#	print "separation: vu= $allvu fgts / total = $nbfgt fgts => @part1 / @part2\n";
#	print "separation: new numbering: @child1 / @child2\n";
#	print "separation: $childping1 / $childping2\n";

};

########################################################################################################
 
