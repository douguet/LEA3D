#!/usr/bin/perl


print"SUBROUTINE SCORE OK \n";

sub score{

	local($file,$nbchildevaluer)=@_;

	@score="";
	foreach $sj (1..$nbchildevaluer){
		$score[$sj]=0;
		$properties[$sj]="";
	};

######################## Fonction TEST
# random score from 0 to 100
#	foreach $sj (1..$nbchildevaluer){
#		$score[$sj]=rand(100);
#		$score[$sj]=sprintf"%3.2f",$score[$sj];
#		print "mol$sj $score[$sj]\n";
#	};
	
########## Function to evaluate

print "COMPUTE ENERGY\n" if($param{VERBOSITY} >= 1);
	
# a molecule must have >= $percentproperties % of the properties to go through docking
$percentproperties=50;
       
if($nbpre_fonction > 0){
		
       	#General properties: separate tmpi.sdf and create tmpi.mol2 for each molecule in sdf
	if($evaluate_prop){
		#not solubility
		if($param{CHARGES} eq 'GAS' || $param{CHARGES} eq 'AM1'){
			chop($tmpprop=`perl $leaexe/properties.pl $file 1 $param{CHARGES} 0`);
			unlink "vector_dipole.pdb" if (-e "vector_dipole.pdb");
		}
		else{	
			chop($tmpprop=`perl $leaexe/properties.pl $file 0 - 0`);
		};

		open(MOL,"<properties.out");
		$sj=1;
       		while(<MOL>){
			@getline=split(' ',$_);
			if($getline[0]!~/^#/){
				$propertiestitle="";
				$propertiesvalue="";
				foreach $k (2..@gettitle-1){
					foreach $prop (0..@fprop-1){
						if($fprop[$prop] eq $gettitle[$k]){
							# composite score score[$sj]
							print "conformer $sj $gettitle[$k] = $getline[$k]\n";
							# give a score between 0 and 1
							$score[$sj]=$score[$sj]+(&composite($getline[$k],$fmin[$prop],$fmax[$prop]))*$fw[$prop];	
							$propertiestitle=$propertiestitle."$gettitle[$k] ";
							$propertiesvalue=$propertiesvalue."$getline[$k] ";
						};
					};
				};
				#To retreive all properties:
				#$properties[$sj]=$properties[$sj]."#".$listtitle.$_;	
				#or to retreive selected properties (in function)
				$properties[$sj]=$properties[$sj]."#".$propertiestitle."\n".$propertiesvalue."\n";
				$sj++;
			}
			else{
				$_ =~s/^#//;
				$listtitle=$_;
				@gettitle=split(' ',$_);
			};	
		};
		close(MOL);
		unlink "properties.out";
	};


	#Chemical Function occurences: use the first molecule of the sdf file (no mol2 file)
	if($evaluate_function){
		$pointeur_function=-1;
		foreach $prop (0..@fprop-1){
			$pointeur_function=$prop if($fprop[$prop] eq "function");
		};
		
		$tmpprop="";
		$tmpprop=&chemfunction($file); # apply at the first molecule only
		#print "tmpprop=$tmpprop\n";
		
		$prop3=0;
		@listfonc=split('_',$listefonction);
		$llistfonc=@listfonc;
		$fract=1/$llistfonc;
		foreach $lfi (0..@listfonc-1){
			if($tmpprop =~ /$listfonc[$lfi]/){
				$prop3=$prop3+$fract;
				print"$listfonc[$lfi] OK ($prop3)\n" if($param{VERBOSITY} >= 1);
			};
		};	
		print "Function (search for $listefonction) Score=$prop3 for all conformers (1 means full match)\n";

		# fill @score for each conformers = same score for each
		foreach $sj (1..$nbchildevaluer){ 
		 	$score[$sj]=$score[$sj]+(&composite($prop3,$fmin[$pointeur_function],$fmax[$pointeur_function]))*$fw[$pointeur_function];
			$properties[$sj]=$properties[$sj]."FUNCTIONS conformer $sj = $tmpprop\n";
		};	
	};

	#Pharmacophore: separate tmpi.sdf and creates tmpi.mol2 for each molecule in sdf 
	if($evaluate_pharm){
		print "Pharmacophore is not implemented in this version\n";
	};

	#120 GC Fingerprint: the conversion by sdfmol2 makes that only the first molecule of $file is analysed
	if($evaluate_finger){
		print "GC Fingerprint is not implemented in this version\n";
	};

	#Shape comparison : keep the best score among all conformers
	if($evaluate_shape){
		print "Shape comparison is not implemented in this version\n";
	};

	#exactmatch using classfgt.pl
	if($evaluate_exactmatch){
		print "Exact match is not implemented in this version\n";
	};
};

	if($nbpre_fonction > 0 && $evaluate_dock){
		print "Docking is not implemented in this version\n";
	}
	elsif($evaluate_dock){
		$flagdock=1;
	};

	if($flagdock){
		print "Docking is not implemented in this version\n";	
	};	

########## Convert total score in percentage

	foreach $sj (1..$nbchildevaluer){
		$score[$sj]=sprintf"%3.2f",($score[$sj]/$sumw)*100;
		#print "$sj $score[$sj]\n";
	};

	@scoretmp="";
	@scoretmp=@score;
	&decreasing_order;
	print "Conformers: decreasing order of scores:\n";
	foreach $sj (1..$nbchildevaluer){
		print "$sj $score[$ranking[$sj]] mol$ranking[$sj]\n";
	};	
	
};

################################################################################
################################################################################

sub decreasing_order{

	@ranking="";
	$debj=0;
	$debj=1 if($scoretmp[0] eq "");
	foreach $sj ($debj..@scoretmp-1){
		$ranking[$sj]=$sj;
	};	
	foreach $sj ($debj..@scoretmp-1){
		$maxclass=$scoretmp[$sj];
		$maxi=$sj;
		foreach $k ($sj..@scoretmp-1){
			if ($maxclass < $scoretmp[$k]){
				$maxclass=$scoretmp[$k];
				$maxi=$k;
			};
		};
		if($maxi != $sj){
			$tmp=$scoretmp[$maxi];
			$scoretmp[$maxi]=$scoretmp[$sj];
			$scoretmp[$sj]=$tmp;
			$tmp=$ranking[$sj];
			$ranking[$sj]=$ranking[$maxi];
			$ranking[$maxi]=$tmp;
		};	
	};
};

################################################################################
#################################################################################

sub composite{

        local($valeur,$compmin,$compmax)=@_;

# Fournit un score entre 0 (unfitted) et 1 (perfect fit)

        $flagcompmin=0;
        $flagcompmax=0;


        if ($compmax eq '' && $compmin ne ''){
                $flagcompmin=1;
        }
        elsif($compmin eq '' && $compmax ne ''){
                $flagcompmax=1;
        };

        $val=0;
        if ($flagcompmin){
                $val=1 if ($valeur >= $compmin);
        }
        elsif($flagcompmax){
                $val=1 if ($valeur <= $compmax);
        }
        else{

		# in order to give a score to value far from required
		# $ecartminmax=($compmax-$compmin)/4; #hard
		$ecartminmax=($compmax-$compmin); #smooth

       	 	$sqecartminmax=$ecartminmax*$ecartminmax;

        	if ($compmin == $compmax){
                	  if ($valeur != $compmin){
	
				# in order to give a score to value far from required
				#   $ecartmin=$compmin*3/4; # hard
				#   0.9 means that val will not be zero from 10% to 100%
				$ecartmin=$compmin*0.9; # smooth give score from 10% of the target value

                        	$sqecartmin=$ecartmin*$ecartmin;
                        	if ($valeur < $compmin){
                                	$val=$sqecartmin-(($compmin-$valeur)*($compmin-$valeur));
                                	if($val < 0){# in case ($compmin-$valeur)**2 is greater than $sqecartmin
						$val=0;
					}
					else{
						$val=$val/$sqecartmin;
					};	
                        	}
                        	else{
                                	$val=$sqecartmin-(($valeur-$compmin)*($valeur-$compmin));
                                	if($val < 0){# in case ($compmin-$valeur)**2 is greater than $sqecartmin
						$val=0;
					}
					else{
						$val=$val/$sqecartmin;
					};	
                        	};
                  	}
                  	else{
                        	$val=1;
                  	};
		}
        	else{
                	if ($valeur < $compmin){
                       		#$val=$sqecartminmax-(($compmin-$valeur)*($compmin-$valeur));
				#$val=$val/$sqecartminmax;
				
                        	$ecartmin=$compmin*0.9;
				$sqecartmin=$ecartmin*$ecartmin;
				$val=$sqecartmin-(($compmin-$valeur)*($compmin-$valeur));

				$val=$val/$sqecartmin;
                  	}
                 	elsif ($valeur > $compmax){
                        	#$val=$sqecartminmax-(($valeur-$compmax)*($valeur-$compmax));
                        	#$val=$val/$sqecartminmax;
				
				$ecartmin=$compmax*0.9;
				$sqecartmin=$ecartmin*$ecartmin;
				$val=$sqecartmin-(($valeur-$compmin)*($valeur-$compmin));

				$val=$val/$sqecartmin;
                  	}
                  	else{# between or == $compmax or == $compmin
                        	$val=1;
                  	};
		};
	};

        if(($val>1)||($val<0)|| $valeur eq ''){
        	$val=0;
        }
	else{
		$val=sprintf"%4.3f",($val);
	};
	$val;
};

#################################################################################
