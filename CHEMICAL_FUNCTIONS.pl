#!/usr/bin/perl


print"SUBROUTINE CHEMICAL_FUNCTIONS \n";

sub chemfunction{

	local($filesdf)=@_;

	# The default van der Waals radii are taken from A. Bondi (1964) "van der Waals Volumes and Radii" J.Phys.Chem. 68, 441-451
	#If an element does not appear in the table it is assigned a value of 2.0?

	#   Ag  1.72      Ar  1.88     As  1.85     Au  1.66
	#   Br  1.85      C   1.70     Cd  1.58     Cl  1.75
	#   Cu  1.40      F   1.47     Ga  1.87     H   1.20
	#   He  1.40      Hg  1.55     I   1.98     In  1.93
	#   K   2.75      Kr  2.02     Li  1.82     Mg  1.73
	#   N   1.55      Na  2.27     Ne  1.54     Ni  1.63
	#   O   1.52      P   1.80     Pb  2.02     Pd  1.63
	#   Pt  1.72      S   1.80     Se  1.90     Si  2.10
	#   Sn  2.17      Te  2.06     Tl  1.96     U  1.86
	#   Xe  2.16      Zn  1.39
	#For boron B used SERF program vdW radii

	%tabR=(
		'C',1.70,
		'O',1.52,
		'N',1.55,
		'S',1.80,
		'P',1.80,
		'B',1.72,
		'Br',1.85,
		'Cl',1.75,
		'I',1.98,
		'F',1.47,
		'H',1.2,
		'Hp',1.1,
	);

        %tabmm=(
                'C',12.01,
                'O',16,
                'N',14.01,
                'S',32.07,
                'P',30.97,
                'B',10.81,
                'Br',79.9,
                'Cl',35.45,
                'I',126.9,
                'F',19,
                'H',1.008,
        );

	##### PROPRIETES FROM SDF
	$stops=0;
	$flagnew=1;
	open(MOL,"<$filesdf");
	while(<MOL>){
	if($stops==0){

		if($flagnew){
			$masse=0;
			$compt=0;
			$ig=1;
			$jg=0;
			@strx='';
			@stry='';
			@strz='';
			@atom='';
			@coval='';
			@fonc='';
			@ifonc='';
			@covfonc='';
			$blanc=' ';	
			@radius='';
			$atomlourd=0;
			$flagnew=0;
			$setchg=0;
		};
		@getstr = split(' ',$_);
		$compt++;
		if (($compt > 4) && ($ig <= $istratom)){
			$strx[$ig]=$getstr[0];
			$stry[$ig]=$getstr[1];
			$strz[$ig]=$getstr[2];
			$atom[$ig]=$getstr[3];
			$atomlourd ++ if($getstr[3] ne 'H');
			$radius[$ig]=$tabR{$getstr[3]};
			$masse=$masse+$tabmm{$getstr[3]};
			$ig++;
		};
		if (($compt > 4) && ($ig > $istratom) && ($jg <=$istrbond)){
			if ($jg == 0){
				$jg++;
			}
			else{
                               @coller=split(' *',$getstr[0]);
                                @coller2=split(' *',$getstr[1]);
                                if(@coller==6 && $getstr[1] ne ""){
                                        $getstr[0]=$coller[0].$coller[1].$coller[2];
                                        $getstr[2]=$getstr[1];
                                        $getstr[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $getstr[1] eq ""){
                                        $getstr[0]=$coller[0].$coller[1];
                                        $getstr[1]=$coller[2].$coller[3].$coller[4];
                                        $getstr[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($_=~/^\s/){
                                                $getstr[0]=$coller[0].$coller[1];
                                                $getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $getstr[0]=$coller[0].$coller[1].$coller[2];
                                                $getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($_=~/^\s/){
                                                $getstr[0]=$coller[0];
                                                $getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $getstr[0]=$coller[0].$coller[1].$coller[2];
						$getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[3];
                                        };
                                }
                                elsif(@coller2==4){
                                        $getstr[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $getstr[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $getstr[0]=$coller[0].$coller[1].$coller[2];
                                        $getstr[1]=$coller[3].$coller[4].$coller[5];
                                        $getstr[2]=$coller[6];
                                };

				$fonc[$getstr[0]]=$fonc[$getstr[0]].$blanc.$getstr[2].'-'.$atom[$getstr[1]].$blanc;
				$ifonc[$getstr[0]]=$ifonc[$getstr[0]].$blanc.$getstr[1].$blanc;
				$covfonc[$getstr[0]]=$covfonc[$getstr[0]].$blanc.$getstr[2];
				$coval[$getstr[0]]=$coval[$getstr[0]]+$getstr[2];

				$fonc[$getstr[1]]=$fonc[$getstr[1]].$blanc.$getstr[2].'-'.$atom[$getstr[0]].$blanc;
				$ifonc[$getstr[1]]=$ifonc[$getstr[1]].$blanc.$getstr[0].$blanc;
				$covfonc[$getstr[1]]=$covfonc[$getstr[1]].$blanc.$getstr[2];
				$coval[$getstr[1]]=$coval[$getstr[1]]+$getstr[2];

				$jg++;
			};
		};
		if ($compt == 4){
			$istratom=$getstr[0];
			$istrbond=$getstr[1];

			@coller=split(' *',$istratom);
			if(@coller>3 && @coller==6){
				$istratom=$coller[0].$coller[1].$coller[2];
				$istrbond=$coller[3].$coller[4].$coller[5];
			}
			elsif(@coller>3 && @coller==5){
				if($_=~/^\s/){
					$istratom=$coller[0].$coller[1];
					$istrbond=$coller[2].$coller[3].$coller[4];
				}
				else{
					$istratom=$coller[0].$coller[1].$coller[2];
					$istrbond=$coller[3].$coller[4];
				};
			};

		};
		if ($_=~/\$\$\$\$/){
			$flagnew=1;
			print "atomes = $istratom\n" if($param{VERBOSITY} >= 1);
			$result_search=&fonction;
			#print "function seen = $result_search\n";
			$stops=1; # read the first molecule only
		};
	};#end if($stops)	
	};
	close(MOL);

	$result_search;
};

########################################################


sub fonction{

## SEARCH FUNCTIONS

	$function='';
	print"$file = " if($param{VERBOSITY} >= 1);

	foreach $l (1..$istratom){
	
	$function=$function.$blanc.'F'.$blanc if($atom[$l] eq 'F');
	$function=$function.$blanc.'Cl'.$blanc if($atom[$l] eq 'Cl');
	$function=$function.$blanc.'Br'.$blanc if($atom[$l] eq 'Br');
	$function=$function.$blanc.'P'.$blanc if($atom[$l] eq 'P');
	$function=$function.$blanc.'I'.$blanc if($atom[$l] eq 'I');

############################################# Carbone

	if($atom[$l] eq 'C'){

		$function=$function.$blanc.'C'.$blanc;

		if($fonc[$l] =~/ 2-O / && $fonc[$l] =~/ 1-O / && $fonc[$l] =~/ 1-N /){
			$function=$function.$blanc.'carbamate'.$blanc;
#			print"carbamate ";
		};
		if($fonc[$l] =~/ 2-O / && $fonc[$l] =~/ 1-O / ){
			$car=$fonc[$l];
			$nocc=$car=~s/ 1-O / 1-O /g;
#			print"nb occ = $nocc\n";
			if ($nocc == 1){ 
				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k (0..@det-1){
					$p=$k if($det[$k] =~/1-O/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				@det=split(' ',$fonc[$iatom]);
				if($fonc[$iatom] =~/ 1-H / || @det==1){
					$function=$function.$blanc.'acid'.$blanc;
#					print"acid ";
				}
				else{
					$function=$function.$blanc.'ester'.$blanc;
#					print"ester ";
				};
			}
			else{
#				print"OC(=O)O ";
			};
		};
		if($fonc[$l] =~/ 2-O / && $fonc[$l] =~/ 1-N / ){
			@det=split(' ',$fonc[$l]);
			$p=0;
			foreach $k (0..@det-1){
				$p=$k if($det[$k] =~/1-N/);
			};
			@det=split(' ',$ifonc[$l]);
			$iatom=$det[$p];
			$car=$fonc[$iatom];
			$nocc=$car=~s/ 1-H / 1-H /g;
			if($nocc == 2){
				$function=$function.$blanc.'amide-ter'.$blanc;
#				print"amide-ter ";
			}
			else{
				$function=$function.$blanc.'amide'.$blanc;
#				print"amide ";
			};
		};

		if($fonc[$l] =~/ 2-O / && $fonc[$l] =~/ 1-H / && $fonc[$l] =~/ 1-C /){
			$function=$function.$blanc.'aldhehyde'.$blanc;
#			print"aldhehyde ";
		};


		if($fonc[$l] =~/ 2-O / && $fonc[$l] =~/ 1-C /){
			$car=$fonc[$l];
			$nocc=$car=~s/ 1-C / 1-C /g;
#			print"nb occ = $nocc\n";
#			print"keto " if ($nocc == 2); 
			$function=$function.$blanc.'keto'.$blanc if ($nocc == 2);
		};

		if($fonc[$l] =~/ 2-O /){
			$function=$function.$blanc.'carbonyl'.$blanc;
		};	

	};

#####################################################################
############################################# Oxygene

	if($atom[$l] eq 'O'){

		$function=$function.$blanc.'O'.$blanc;

		if($fonc[$l] =~/ 1-C /){
			$car=$fonc[$l];
			$nocc=$car=~s/ 1-C / 1-C /g;
			if($nocc == 2){
				$ok=1;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k (0..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O / || $fonc[$iatom] =~/ 2-S /);
				$pold=$p;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k ($pold+1..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O /  || $fonc[$iatom] =~/ 2-S /);

				if($ok){
					$function=$function.$blanc.'ether'.$blanc;
#					print "ether ";
				};
			}
			elsif($nocc == 1 && $fonc[$l] =~/ 1-H /){	

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k (0..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$nocc2=$car=~s/ 1-H / 1-H /g;
                                if($nocc2 == 2){
                                        $function=$function.$blanc.'alcohol1'.$blanc;
#                                       print"alcohol1 ";
                                }
                                elsif($nocc2 == 1 && $fonc[$iatom] !~/ 2-O / && $fonc[$iatom] !~/ 2-S /){
                                        $function=$function.$blanc.'alcohol2'.$blanc;
#                                       print"alcohol2 ";
                                }
                                elsif($nocc2 == 0 && $fonc[$iatom] !~/ 2-O / && $fonc[$iatom] !~/ 2-S /){
                                        $function=$function.$blanc.'alcohol3'.$blanc;
#                                       print"alcohol3 ";
                                };
			};
		};
	};

#####################################################################
############################################# S

	if($atom[$l] eq 'S'){
		$function=$function.$blanc.'S'.$blanc;
		if($fonc[$l] =~/ 1-C /){
			$car=$fonc[$l];
			$nocc=$car=~s/ 1-C / 1-C /g;
			if($nocc == 1 && $fonc[$l] =~/ 1-H / && $coval[$l] == 2){
				$function=$function.$blanc.'thiol'.$blanc;
			};
		};		 
	};

#####################################################################
############################################# Azote

	if($atom[$l] eq 'N'){
		$function=$function.$blanc.'N'.$blanc;
		if($fonc[$l] =~/ 1-C /){
			$car=$fonc[$l];
			$nocc=$car=~s/ 1-C / 1-C /g;
			if($nocc == 3){
				$ok=1;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k (0..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O / || $fonc[$iatom] =~/ 2-S /);
				$pold=$p;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k ($pold+1..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O / || $fonc[$iatom] =~/ 2-S /);
				$pold=$p;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k ($pold+1..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O/ || $fonc[$iatom] =~/ 2-S /);

				if($ok){
					$function=$function.$blanc.'amine3'.$blanc;
#					print"amine3 ";
				};
				
			}
			elsif($nocc == 2 && $fonc[$l] =~/ 1-H /){
				$ok=1;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k (0..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O / || $fonc[$iatom] =~/ 2-S /);
				$pold=$p;

				@det=split(' ',$fonc[$l]);
				$p=0;
				foreach $k ($pold+1..@det-1){
					$p=$k if($det[$k] =~/1-C/);
				};
				@det=split(' ',$ifonc[$l]);
				$iatom=$det[$p];
				$car=$fonc[$iatom];
				$ok=0 if($fonc[$iatom] =~/ 2-O / || $fonc[$iatom] =~/ 2-S /);

				if($ok){
					$function=$function.$blanc.'amine2'.$blanc;
#					print"amine2 ";
				};
			}
			elsif($nocc == 1 && $fonc[$l] =~/ 1-H /){	
				$car=$fonc[$l];
				$nocc2=$car=~s/ 1-H / 1-H /g;
				if($nocc2 == 2){
					$ok=1;

					@det=split(' ',$fonc[$l]);
					$p=0;
					foreach $k (0..@det-1){
						$p=$k if($det[$k] =~/1-C/);
					};
					@det=split(' ',$ifonc[$l]);
					$iatom=$det[$p];
					$car=$fonc[$iatom];
					$ok=0 if($fonc[$iatom] =~/ 2-O / || $fonc[$iatom] =~/ 2-S /);

					if($ok){
						$function=$function.$blanc.'amine1'.$blanc;
#						print"amine1 ";
					};
				};
			};
		};
	};
#####################################################################

	};


	$prop2='';

	$car=$function;
	$nocc=$car=~s/ C / C /g;
	$prop2=$prop2.'_'.$nocc.'C' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ N / N /g;
        $prop2=$prop2.'_'.$nocc.'N' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ O / O /g;
        $prop2=$prop2.'_'.$nocc.'O' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ S / S /g;
        $prop2=$prop2.'_'.$nocc.'S' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ P / P /g;
        $prop2=$prop2.'_'.$nocc.'P' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ F / F /g;
        $prop2=$prop2.'_'.$nocc.'F' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ Br / Br /g;
        $prop2=$prop2.'_'.$nocc.'Br' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ Cl / Cl /g;
        $prop2=$prop2.'_'.$nocc.'Cl' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ I / I /g;
        $prop2=$prop2.'_'.$nocc.'I' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ acid / acid /g;
	$prop2=$prop2.'_'.$nocc.'acid' if ($nocc > 0);
	
	$car=$function;
	$nocc=$car=~s/ ester / ester /g;
	$prop2=$prop2.'_'.$nocc.'ester' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ carbonyl / carbonyl /g;
	$prop2=$prop2.'_'.$nocc.'carbonyl' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ carbamate / carbamate /g;
	$prop2=$prop2.'_'.$nocc.'carbamate' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ amide / amide /g;
	$prop2=$prop2.'_'.$nocc.'amide' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ amide-ter / amide-ter /g;
	$prop2=$prop2.'_'.$nocc.'amide-ter' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ aldhehyde / aldhehyde /g;
	$prop2=$prop2.'_'.$nocc.'aldhehyde' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ keto / keto /g;
	$prop2=$prop2.'_'.$nocc.'keto' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ amine1 / amine1 /g;
	$prop2=$prop2.'_'.$nocc.'amine1' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ amine2 / amine2 /g;
	$prop2=$prop2.'_'.$nocc.'amine2' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ amine3 / amine3 /g;
	$prop2=$prop2.'_'.$nocc.'amine3' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ alcohol1 / alcohol1 /g;
	$prop2=$prop2.'_'.$nocc.'alcohol1' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ alcohol2 / alcohol2 /g;
	$prop2=$prop2.'_'.$nocc.'alcohol2' if ($nocc > 0);

        $car=$function;
        $nocc=$car=~s/ alcohol3 / alcohol3 /g;
        $prop2=$prop2.'_'.$nocc.'alcohol3' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ thiol / thiol /g;
	$prop2=$prop2.'_'.$nocc.'thiol' if ($nocc > 0);

	$car=$function;
	$nocc=$car=~s/ ether / ether /g;
	$prop2=$prop2.'_'.$nocc.'ether' if ($nocc > 0);

	print "FUNCTIONS: $prop2\n";# if($param{VERBOSITY} >= 1);
	$prop2;

};


#################################################################################################"
#################################################################################################"

