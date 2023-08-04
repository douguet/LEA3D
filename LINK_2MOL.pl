#!/usr/bin/perl

### Combine deux fragments et ecrit la structure dans file_out
### link2mol(file1,atom1,file2,atom2,option) 
#   Move file1 only to connect file2 if $origin==0 or ==2 else don't move (option 1 that uses <ORIGIN>]
#   $origin=0 use number of atoms to link fragments or use the first X if number of atoms are set to 0
#   $origin=1 use POINTS and ORIGIN datablocks and do no move fgts (number of atoms are set to 0)
#   $origin=2 in cases of lego (LEA3D) ; print datablock LEGO 
#   New molecule numbering begins with atoms from file1 and then atoms from file2

	local($file1,$noatom1,$file2,$noatom2,$origin)=@ARGV;

	if($file1 eq '' || $file2 eq '' || $noatom1 eq '' || $noatom2 eq ''){
		die "usage: link2mol <file1.sdf> <no_atom_1 or '0'> <file2.sdf> <no_atom_2 or '0'>  <(optional) set '1' to use ORIGIN datablock to re-build (fgts do not move)>\nfile2.sdf is fixed during the process\n"; 
	};

	$origin=0 if($origin eq "");
	# if $origin==2 then launched by eDesign to connect lego
	 
	$file3 = "combin.sdf";
	$moli=1;

###################################################### REPRIS DE LINK

## RECHERCHE FGT 1

# MODIF
	$no1=1;
##### FIN MODIF

@lego1="";
$lego1i=0;
$pointsfgt1='';
$flagnew=1;
$no=0;
open(IN,"<$file1");
while (<IN>){
	$conv2=$_;
	if ($flagnew){
		$compt=0;
		$i=1;
		$j=0;
		$blanc=" ";
		$flagnew=0;
		$no++;
		if($no == $no1){
			@type1='';
			@corx1='';
			@cory1='';
			@corz1='';
			@lignecor1='';
			@lignebond1='';
			$lignedata='';
			@fonc1='';
			@ifonc1='';
			@covfonc1='';
			$posl=0;
			$fposl=0;
			$posr=0;
			$fposr=0;
			$regardfgt=0;
			$regardlego=0;
			$flagorigin=0;
		};
	};
	if($no == $no1){
		@get = split(' ',$_);
		$compt++;
	
		if($flagorigin){
			$lignedata=$get[0];
			$flagorigin=0;
		};
			
		if ($fposl){
			$posl=$get[0];
			$fposl=0;
		};
		if ($fposr){
			$posr=$get[0];
			$fposr=0;
		};

		if (($compt > 4) && ($i <= $nbatom1)){
			$corx1[$i]=$get[0];
			$cory1[$i]=$get[1];
			$corz1[$i]=$get[2];
			$type1[$i]=$get[3];
			$lignecor1[$i]=$_;
			$i++;
		};
		if (($compt > 4) && ($i > $nbatom1) && ($j <= $nbond1)){
			if ($j == 0){
				$j++;
			}
			else{
                               @coller=split(' *',$get[0]);
                                @coller2=split(' *',$get[1]);
                                if(@coller==6 && $get[1] ne ""){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[2]=$get[1];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $get[1] eq ""){
                                        $get[0]=$coller[0].$coller[1];
                                        $get[1]=$coller[2].$coller[3].$coller[4];
                                        $get[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($_=~/^\s/){
                                                $get[0]=$coller[0].$coller[1];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($_=~/^\s/){
                                                $get[0]=$coller[0];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
						$get[2]=$get[1];
                                                $get[1]=$coller[3];
                                        };					
                                }
                                elsif(@coller2==4){
                                        $get[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $get[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                        $get[2]=$coller[6];
                                };

		$fonc1[$get[0]]=$fonc1[$get[0]].$blanc.$get[2].'-'.$type1[$get[1]].$blanc;
		$ifonc1[$get[0]]=$ifonc1[$get[0]].$blanc.$get[1].$blanc;
		$covfonc1[$get[0]]=$covfonc1[$get[0]].$blanc.$get[2] if($type1[$get[1]] ne 'H' && $type1[$get[1]] ne 'X');


		$fonc1[$get[1]]=$fonc1[$get[1]].$blanc.$get[2].'-'.$type1[$get[0]].$blanc;
		$ifonc1[$get[1]]=$ifonc1[$get[1]].$blanc.$get[0].$blanc;
		$covfonc1[$get[1]]=$covfonc1[$get[1]].$blanc.$get[2] if($type1[$get[0]] ne 'H' && $type1[$get[0]] ne 'X');

			$lignebond1[$j]=$_;
			$j++;
			};
		};
		if ($compt == 4){
			$nbatom1=$get[0];
			$nbond1=$get[1];

			@coller=split(' *',$nbatom1);
			if(@coller > 3 && @coller == 6){
				$nbatom1=$coller[0].$coller[1].$coller[2];
				$nbond1=$coller[3].$coller[4].$coller[5];
			}
			elsif( @coller > 3 && @coller == 5){
				if($_=~/^\s/){
					$nbatom1=$coller[0].$coller[1];
					$nbond1=$coller[2].$coller[3].$coller[4];
				}
				else{
					$nbatom1=$coller[0].$coller[1].$coller[2];
					$nbond1=$coller[3].$coller[4];
				};
			};
			#print "1_ $nbatom1, $nbond1\n" if($param{VERBOSITY} >= 2);
			#$compt++;
		};

		if ($j > $nbond1 && $_=~/^>/ && $_=~/ORIGIN/){
			$flagorigin=1;
		}; 		
                if($regardfgt){
	        	$pointsfgt1=$get[0];
	        	$regardfgt=0;
	  	};
	        $regardfgt=1 if($_=~/^>/ && $_=~/POINTS/);

		$regardlego=0 if($get[0] eq "" || $_=~/^>/);
		if($regardlego){
			$lego1[$lego1i]="$get[0] $get[1] $get[2]";
			$lego1i++;
		};	
		$regardlego=1 if($_=~/^>/ && $_=~/LEGO/);
		
		$fposl=1 if ($_=~/^>/ && $_=~/TAG_l/);
		$fposr=1 if ($_=~/^>/ && $_=~/TAG_r/);

	};
	if ($conv2 =~/\$\$\$\$/){
		$flagnew=1;
	};
};
close(IN);

#print "lego1 @lego1\n";


######################################
## RECHERCHE FGT 2

# MODIF
	$no2=1;
##### FIN MODIF


@lego2="";
$lego2i=0;
$pointsfgt2='';
$flagnew=1;
$no=0;
open(IN,"<$file2");
while (<IN>){
	$conv2=$_;
	if ($flagnew){
		$compt=0;
		$i=1;
		$j=0;
		$blanc=" ";
		$flagnew=0;
		$no++;
		if($no == $no2){
			@type2='';
			@corx2='';
			@cory2='';
			@corz2='';
			@lignecor2='';
			@lignebond2='';
			$lignedata2="";
			@fonc2='';
			@ifonc2='';
			@covfonc2='';
			$posl='';
			$fposl=0;
			$posr='';
			$fposr=0;
			$regardlego=0;
			$regardfgt=0;
			$flagorigin=0;
		};

	};
	if($no == $no2){
		@get = split(' ',$_);
		$compt++;

                if($flagorigin){
                        $lignedata2=$get[0];
			$flagorigin=0;
                };

		if ($fposl){
			$posl=$get[0];
			$fposl=0;
		};
		if ($fposr){
			$posr=$get[0];
			$fposr=0;
		};
		

		if (($compt > 4) && ($i <= $nbatom2)){
			$corx2[$i]=$get[0];
			$cory2[$i]=$get[1];
			$corz2[$i]=$get[2];
			$type2[$i]=$get[3];
			$lignecor2[$i]=$_;
			$i++;
		};
		if (($compt > 4) && ($i > $nbatom2) && ($j <= $nbond2)){
			if ($j == 0){
				$j++;
			}
			else{

                               @coller=split(' *',$get[0]);
                                @coller2=split(' *',$get[1]);
                                if(@coller==6 && $get[1] ne ""){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[2]=$get[1];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $get[1] eq ""){
                                        $get[0]=$coller[0].$coller[1];
                                        $get[1]=$coller[2].$coller[3].$coller[4];
                                        $get[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($_=~/^\s/){
                                                $get[0]=$coller[0].$coller[1];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($_=~/^\s/){
                                                $get[0]=$coller[0];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
						$get[2]=$get[1];
                                                $get[1]=$coller[3];
                                        };					
                                }
                                elsif(@coller2==4){
                                        $get[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $get[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                        $get[2]=$coller[6];
                                };

		$fonc2[$get[0]]=$fonc2[$get[0]].$blanc.$get[2].'-'.$type2[$get[1]].$blanc;
		$ifonc2[$get[0]]=$ifonc2[$get[0]].$blanc.$get[1].$blanc;
		$covfonc2[$get[0]]=$covfonc2[$get[0]].$blanc.$get[2] if($type2[$get[1]] ne 'H' && $type2[$get[1]] ne 'X');
	

		$fonc2[$get[1]]=$fonc2[$get[1]].$blanc.$get[2].'-'.$type2[$get[0]].$blanc;
		$ifonc2[$get[1]]=$ifonc2[$get[1]].$blanc.$get[0].$blanc;
		$covfonc2[$get[1]]=$covfonc2[$get[1]].$blanc.$get[2] if($type2[$get[0]] ne 'H' && $type2[$get[0]] ne 'X');

			$lignebond2[$j]=$_;
			$j++;
			};
		};

		if ($compt == 4){
			$nbatom2=$get[0];
			$nbond2=$get[1];
		
			@coller=split(' *',$nbatom2);
			if(@coller>3 && @coller==6){
				$nbatom2=$coller[0].$coller[1].$coller[2];
				$nbond2=$coller[3].$coller[4].$coller[5];
			}
			elsif(@coller>3 && @coller==5){
				if($_=~/^\s/){
					$nbatom2=$coller[0].$coller[1];
					$nbond2=$coller[2].$coller[3].$coller[4];
				}
				else{
					$nbatom2=$coller[0].$coller[1].$coller[2];
					$nbond2=$coller[3].$coller[4];
				};
			};
			#print "2_ $nbatom2, $nbond2\n" if($param{VERBOSITY} >= 2);
			$compt++;
		};
	
                if ($j > $nbond2 && $_=~/^>/ && $_=~/ORIGIN/){
                        $flagorigin=1;
                };

	 	if($regardfgt){
	        	$pointsfgt2=$get[0];
	        	$regardfgt=0;
	  	};
	        $regardfgt=1 if($_=~/^>/ && $_=~/POINTS/);
	
		$regardlego=0 if($get[0] eq "" || $_=~/^>/);
		if($regardlego){
			$lego2[$lego2i]="$get[0] $get[1] $get[2]";
			$lego2i++;
		};
		$regardlego=1 if($_=~/^>/ && $_=~/LEGO/);
															
		$fposl=1 if ($_=~/^>/ && $_=~/TAG_l/);
		$fposr=1 if ($_=~/^>/ && $_=~/TAG_r/);
	};

	if ($conv2 =~/\$\$\$\$/){
		$flagnew=1;
	};
};
close(IN);

#print "lego2 @lego2\n";
#if($origin){
	#print "ORIGIN $lignedata and POINTS $pointsfgt1\n";
	#print "ORIGIN $lignedata2 and POINTS $pointsfgt2\n";
#};

###########################################################

$l1="";
$l2="";
#if $pointsfgt1 ne ""
$l1i="";
$l2i="";

if($origin==1){

	@spptfgto1=split('-',$lignedata);
	@spptfgt1=split('-',$pointsfgt1);
	@spptfgto2=split('-',$lignedata2);
	@spptfgt2=split('-',$pointsfgt2);

	$gorigin=0;
	foreach $bo (0..@spptfgto1-1){
		@spptfgto12=split('_',$spptfgto1[$bo]);
		foreach $bo2 (0..@spptfgto2-1){	
			@spptfgto22=split('_',$spptfgto2[$bo2]);
			if(($spptfgto22[0] == $spptfgto12[0] && $spptfgto22[1] == $spptfgto12[1]) || ($spptfgto22[0] == $spptfgto12[1] && $spptfgto22[1] == $spptfgto12[0])){
				$gorigin=1;
				#print "$bo $spptfgto1[$bo] $bo2 $spptfgto2[$bo2] idem\n";
				$point1=$bo;
				$point2=$bo2;
				last;
			};
		};
		last if($gorigin);
	};
	if($gorigin){
		#print "idem $point1 $spptfgto1[$point1] et $point2 $spptfgto2[$point2]\n";		
		$l1=$spptfgt1[$point1];
		$l2=$spptfgt2[$point2];
		$l1i=$point1;
		$l2i=$point2;
	};

}
else{
	if($noatom1==0 && $noatom2==0 && $pointsfgt1 ne "" && $pointsfgt2 ne ""){
		@spptfgt1=split('-',$pointsfgt1);
		@spptfgt2=split('-',$pointsfgt2);
		$l1=$spptfgt1[0];
		$l2=$spptfgt2[0];
		$l1i=0;
		$l2i=0;
	}
	else{
		$l1=$noatom1;
		$l2=$noatom2;
		$l1i=-1;#LEA3D - impossible no of atom
		$l2i=-1;
	};	
};

if($l1 eq "" || $l2 eq ""){
	print "no atoms (l1=$l1 and l2=$l2) to link\n";
	exit(0);
};

        @det=split(' ',$fonc1[$l1]);
        $p="";
        $px="";
        @tabiatom1="";
        $tabia1=0;
        foreach $k (0..@det-1){
                $p=$k if($det[$k] =~/1-H/);
                $px=$k if($det[$k] =~/1-X/);
                if($det[$k] =~/1-X/){
                        @det2=split(' ',$ifonc1[$l1]);
                        $tabiatom1[$tabia1]=$det2[$k];
                        $tabia1++;
                };
        };
        @det=split(' ',$ifonc1[$l1]);
        $iatom1=$det[$p];
        $iatom1x=$det[$px];
        $iatom1=$iatom1x if($px ne "");

	#print "$type1[$l1] $l1 -- $type1[$iatom1] $iatom1 sera substitue\n" ;

	@det=split(' ',$fonc2[$l2]);
	$p="";
	$px="";
	#$p2=0;
	@tabiatom2="";
	$tabia2=0;
	foreach $k (0..@det-1){
		$p=$k if($det[$k] =~/1-H/);
		$px=$k if($det[$k] =~/1-X/);
		#$p2=$k if($det[$k] =~/1-X/);
                if($det[$k] =~/1-X/){
			@det2=split(' ',$ifonc2[$l2]);
                        $tabiatom2[$tabia2]=$det2[$k];
                        $tabia2++;
                };

	};
	@det=split(' ',$ifonc2[$l2]);
	$iatom2=$det[$p];
	$iatom2x=$det[$px];
	#print "$det[$px]\n";
	#$iatom2=$det[$p2];
	$iatom2=$iatom2x if($px ne "");
	#print "$iatom2\n";
	
	#print "$type2[$l2] $l2 -- $type2[$iatom2] $iatom2 sera substitue\n";

#################################################################################
# Alignement des vecteurs
################################################################################
if($origin ==1 && $tabiatom1[0] ne "" && $tabiatom2[0] ne ""){

	foreach $tabia1 (0..@tabiatom1-1){
		foreach $tabia2 (0..@tabiatom2-1){

			#print "$tabia1 ($tabiatom1[$tabia1]) $tabia2 ($tabiatom2[$tabia2])\n";
 
			$iatom1=$tabiatom1[$tabia1];
			$iatom2=$tabiatom2[$tabia2];
	
########################################
# Vecteurs

#print"corx2[iatom2] = $corx2[$iatom2],$cory2[$iatom2],$corz2[$iatom2]\n";

#print"$corx1[$iatom1] (= 0 en principe)\n";
	$ux=$corx1[$l1]-$corx1[$iatom1];
	$uy=$cory1[$l1]-$cory1[$iatom1];
	$uz=$corz1[$l1]-$corz1[$iatom1];
	$nu=sqrt($ux**2+$uy**2+$uz**2);
#print"nu = $nu\n";

#print"$corx2[$l2] (= 0 en principe)\n";
	$vx=$corx2[$iatom2]-$corx2[$l2];
	$vy=$cory2[$iatom2]-$cory2[$l2];
	$vz=$corz2[$iatom2]-$corz2[$l2];
	$nv=sqrt($vx**2+$vy**2+$vz**2);
#print"nv = $nv\n";

########################################
# Angle entre les deux vecteurs a fusionner


	$cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);

	$testcos=1-$cosangle*$cosangle;
#print"$testcos\n";
	$testcos=sqrt($testcos*$testcos);
#print"$testcos\n";
	$teta=atan2(sqrt($testcos),$cosangle);
	$tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;

#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz)\n";

			if($tetab <  0.5){
				last;
			};

		};
		if($tetab <  0.5){
			last;
		};
	};
	
	if($tetab > 0.5){
	#if($tetab > 30){#modify may 2020 to build bioisoster replacement
	#die "in ORIGIN option, Angle(vecteur1,vecteur2) $tetab doit etre proche de 0 !!\n";
		print "warning: in ORIGIN option, Angle(vecteur1,vecteur2) $tetab must be close to 0.0\n";
	};

}
else{

if($origin==0 || $origin==2){

	## centre de masse en iatom1
        $corx=$corx1[$iatom1];
        $cory=$cory1[$iatom1];
        $corz=$corz1[$iatom1];
        foreach $k (1..$nbatom1){
                $corx1[$k]=$corx1[$k]-$corx;
                $cory1[$k]=$cory1[$k]-$cory;
                $corz1[$k]=$corz1[$k]-$corz;
        };

	########  TRANSLATION sur iatom1

	#print "H sub = $corx1[$iatom1],$cory1[$iatom1],$corz1[$iatom1]\n";
	#print "1_C = $corx2[$l2],$cory2[$l2],$corz2[$l2]\n";
        $dx=$corx2[$l2]-$corx1[$iatom1];
        $dy=$cory2[$l2]-$cory1[$iatom1];
        $dz=$corz2[$l2]-$corz1[$iatom1];

	#print "$dx,$dy,$dz\n";

        foreach $j (1..$nbatom2){
                $corx2[$j]=$corx2[$j]-$dx;
                $cory2[$j]=$cory2[$j]-$dy;
                $corz2[$j]=$corz2[$j]-$dz;
        };

	#print "2_C = $corx2[$l2],$cory2[$l2],$corz2[$l2]\n";

	#print"translation $dx, $dy, $dz\n";

#################################################################################
# Alignement des vecteurs
################################################################################

########################################
# Vecteurs

	#print"corx2[iatom2] = $corx2[$iatom2],$cory2[$iatom2],$corz2[$iatom2]\n";

	#print"$corx1[$iatom1] (= 0 en principe)\n";
        $ux=$corx1[$l1]-$corx1[$iatom1];
        $uy=$cory1[$l1]-$cory1[$iatom1];
        $uz=$corz1[$l1]-$corz1[$iatom1];
        $nu=sqrt($ux**2+$uy**2+$uz**2);
	#print"nu = $nu\n";

	#print"$corx2[$l2] (= 0 en principe)\n";
        $vx=$corx2[$iatom2]-$corx2[$l2];
        $vy=$cory2[$iatom2]-$cory2[$l2];
        $vz=$corz2[$iatom2]-$corz2[$l2];
        $nv=sqrt($vx**2+$vy**2+$vz**2);
	#print"nv = $nv\n";

########################################
# Angle entre les deux vecteurs a fusionner u (file1) and v (file2)

        $cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);

        $testcos=1-$cosangle*$cosangle;
	#print"$testcos\n";
        $testcos=sqrt($testcos*$testcos);
	#print"$testcos\n";
        $teta=atan2(sqrt($testcos),$cosangle);
        $tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;

	#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz)\n";
	#print"AVANT angle(vecteur1,vecteur2) $tetab\n";

########################################
# Vecteur normes

        $ux=$ux/$nu;
        $uy=$uy/$nu;
        $uz=$uz/$nu;
	#print"$ux,$uy,$uz\n";

        $vx=$vx/$nv;
        $vy=$vy/$nv;
        $vz=$vz/$nv;
	#print"$vx,$vy,$vz\n";

########################################
# Calcul des vecteurs orthogonaux

# W :
        $wx=($uy*$vz)-($uz*$vy);
        $wy=($uz*$vx)-($ux*$vz);
        $wz=($ux*$vy)-($uy*$vx);
	#print"$wx,$wy,$wz\n";
        $nw=sqrt($wx**2+$wy**2+$wz**2);
	#print"nw = $nw\n";

	#IF VECTORS ARE NOT ALIGNED:
	if($nw != 0){

        	$wx=$wx/$nw;
        	$wy=$wy/$nw;
        	$wz=$wz/$nw;

	# K :
	        $kx=($wy*$uz)-($wz*$uy);
        	$ky=($wz*$ux)-($wx*$uz);
        	$kz=($wx*$uy)-($wy*$ux);
        	$nk=sqrt($kx**2+$ky**2+$kz**2);
		#print"nk = $nk\n";
        	$kx=$kx/$nk;
        	$ky=$ky/$nk;
        	$kz=$kz/$nk;

########################################
# ROTATION
		# valeur de l'angle between v (file2) and k (orthogonal vector between u (file1) and w (orthogonal u,v))

        	$cosangle=(($vx*$kx)+($vy*$ky)+($vz*$kz))/($nv * $nk);
       		
		#$angle=atan2(sqrt(1-$cosangle*$cosangle),$cosangle)*180/3.14159;
       		#print"angle(v,k) $angle (cos=$cosangle)\n";
        
		#??which teta angle ?? the one from vector u and v ??
		$teta=sqrt($teta*$teta);
        	$teta=-1*$teta if($cosangle < 0);

		#print"Rotation angle $tetab\n";
		
		#changement de repere

        	$x=$ux * $ux + $uy * $uy + $uz * $uz;
        	$y=$ux * $kx + $uy * $ky + $uz * $kz;
        	$z=$ux * $wx + $uy * $wy + $uz * $wz;

		# Rotation

        	$xrot= $x*cos($teta) - $y*sin($teta);
        	$yrot= $x*sin($teta) + $y*cos($teta);
        	$zrot=$z;

		# rechangement de repere

        	$x1=$xrot*$ux + $yrot*$kx + $zrot*$wx;
        	$y1=$xrot*$uy + $yrot*$ky + $zrot*$wy;
        	$z1=$xrot*$uz + $yrot*$kz + $zrot*$wz;

#########################################
#24Jan2017 comments:	#only atom defining vector u has to be *$nu because $x, $y and $z have been calculated with $du (normalized instead of $corx1[$l1] .... 
 #       $corx1[$l1]=$x1*$nu;
 #       $cory1[$l1]=$y1*$nu;
 #       $corz1[$l1]=$z1*$nu;
#
#       print"corx1[l1] $corx1[$l1],$cory1[$l1],$corz1[$l1]\n";
#$cosangle=(($corx1[$l1]*$corx2[$iatom2])+($cory1[$l1]*$cory2[$iatom2])+($corz1[$l1]*$corz2[$iatom2]))/(sqrt($corx1[$l1]**2+$cory1[$l1]**2+$corz1[$l1]**2))/(sqrt($corx2[$iatom2]**2+$cory2[$iatom2]**2+$corz2[$iatom2]**2));
#$angle=atan2(sqrt(1-$cosangle*$cosangle),$cosangle)*180.0/3.14159;
#print"vecteur1($corx1[$l1],$cory1[$l1],$corz1[$l1]),vecteur2($corx2[$iatom2],$cory2[$iatom2],$corz2[$iatom2])\n";
#print"APRES angle(vecteur1,vecteur2) $angle\n";
#


#######################################
# pour toutes les autres coordonnees :
# On a besoin de teta, u, k et w

		#print "teta = $teta\n";
		#print "u = $ux,$uy,$uz\n";
		#print "k = $kx,$ky,$kz\n";
		#print "w = $wx,$wy,$wz\n";

		foreach $o (1..$nbatom1){

			#print"$corx1[$o],$cory1[$o],$corz1[$o]\n";

        		#if($o != $l1 && $o != $iatom1){
			if($o != $iatom1){#$iatom1 in (0,0,0)

                		$ncorx=$corx1[$o];
                		$ncory=$cory1[$o];
                		$ncorz=$corz1[$o];

				# changement de repere
                		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
                		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
                		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;

				# Rotation
                		$xrot= $x*cos($teta) - $y*sin($teta);
                		$yrot= $x*sin($teta) + $y*cos($teta);
		                $zrot=$z;

				# rechangement de repere
                		$corx1[$o]=($xrot*$ux + $yrot*$kx + $zrot*$wx);
                		$cory1[$o]=($xrot*$uy + $yrot*$ky + $zrot*$wy);
                		$corz1[$o]=($xrot*$uz + $yrot*$kz + $zrot*$wz);
        		};
			#print"$corx1[$o],$cory1[$o],$corz1[$o]\n";
			#print"\n";
		}; 
	}
	else{
	#IF VECTORS ARE ALIGNED:

        	foreach $o (1..$nbatom1){
		# inverse les coordonnees
		#Nothing to do perfect position : comments 24Jan2017
                	#$corx1[$o]=$corx1[$o]*-1;
                	#$cory1[$o]=$cory1[$o]*-1;
                	#$corz1[$o]=$corz1[$o]*-1;
        	};
	};

	foreach $o (1..$nbatom1){
		$corx1[$o]=$corx1[$o]+$dx;
		$cory1[$o]=$cory1[$o]+$dy;
		$corz1[$o]=$corz1[$o]+$dz;
	};

	foreach $o (1..$nbatom2){
		$corx2[$o]=$corx2[$o]+$dx;
		$cory2[$o]=$cory2[$o]+$dy;
		$corz2[$o]=$corz2[$o]+$dz;
	};

}
else{ #$origin != 0 and $origin != 2
       ## Never entered ??

########################################
# Vecteurs

#print"corx2[iatom2] = $corx2[$iatom2],$cory2[$iatom2],$corz2[$iatom2]\n";

#print"$corx1[$iatom1] (= 0 en principe)\n";
        $ux=$corx1[$l1]-$corx1[$iatom1];
        $uy=$cory1[$l1]-$cory1[$iatom1];
        $uz=$corz1[$l1]-$corz1[$iatom1];
        $nu=sqrt($ux**2+$uy**2+$uz**2);
#print"nu = $nu\n";

#print"$corx2[$l2] (= 0 en principe)\n";
        $vx=$corx2[$iatom2]-$corx2[$l2];
        $vy=$cory2[$iatom2]-$cory2[$l2];
        $vz=$corz2[$iatom2]-$corz2[$l2];
        $nv=sqrt($vx**2+$vy**2+$vz**2);
#print"nv = $nv\n";

########################################
# Angle entre les deux vecteurs a fusionner

        $cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);

        $testcos=1-$cosangle*$cosangle;
	#print"$testcos\n";
        $testcos=sqrt($testcos*$testcos);
	#print"$testcos\n";
        $teta=atan2(sqrt($testcos),$cosangle);
        $tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;

	#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz) angle= $tetab\n";
	if($tetab > 0.5){
        	die "Angle(vecteur1,vecteur2) $tetab doit etre proche de 0 !!\n";
	};

}; # end of if origin==0 or origin==2

};# end if($origin ==1 && $tabiatom1[0] ne "" && $tabiatom2[0] ne ""){


#########################################
# PRINTING is common to all options
# 
	$nbatom3=$nbatom1 + $nbatom2 - 2;
	$nbond3=$nbond1 + $nbond2 -1;

## Print HEAD

	open(OUTC,">$file3");
	printf OUTC "\n";
	printf OUTC "lea \n";
	printf OUTC "\n";
	printf OUTC "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$nbatom3,$nbond3; 	

## Print COORDONNEES
	
	#2021: rdkit is picky with N+ atoms : must add "M  CHG  1  20   1"
	@rdkatom="";
	@rdkvalence="";
	$rdki=1;

	foreach $l (1..$nbatom1){
		@get=split(' ',$lignecor1[$l]);
		if($l != $iatom1){
			$corx1[$l]=sprintf "%4.4f",$corx1[$l];
			$cory1[$l]=sprintf "%4.4f",$cory1[$l];
			$corz1[$l]=sprintf "%4.4f",$corz1[$l];

			#$get[3]=~s/X/H/;   # POUR FGTS avec X au lieu de H

			#set to 0 because Corina uses $get[6] to set stereo instead of 3D coordinates
                        $get[6]=0;

			@tget3=split(' *',$get[3]);
			$lget3=@tget3;
			if($lget3 == 1){
				printf OUTC "%10s%10s%10s%1s%1s%1s%3s%3s%3s%3s%3s%3s\n",$corx1[$l],$cory1[$l],$corz1[$l],$blanc,$get[3],$blanc,$get[4],$get[5],$get[6],$get[4],$get[5],$get[6];
			}
			else{
				printf OUTC "%10s%10s%10s%1s%2s%3s%3s%3s%3s%3s%3s\n",$corx1[$l],$cory1[$l],$corz1[$l],$blanc,$get[3],$get[4],$get[5],$get[6],$get[4],$get[5],$get[6];
			};
			$rdkatom[$rdki]=$get[3];
			$rdki++;
		};
	 }; 

	foreach $l (1..$nbatom2){ 	
		@get=split(' ',$lignecor2[$l]); 
		if($l != $iatom2){
			$corx2[$l]=sprintf "%4.4f",$corx2[$l];
			$cory2[$l]=sprintf "%4.4f",$cory2[$l];
			$corz2[$l]=sprintf "%4.4f",$corz2[$l];

			#$get[3]=~s/X/H/;   # POUR FGTS avec X au lieu de H
			
			#set to 0 because Corina uses $get[6] to set stereo instead of 3D coordinates
			$get[6]=0;

			@tget3=split(' *',$get[3]);
			$lget3=@tget3;
			if($lget3 == 1){
				printf OUTC "%10s%10s%10s%1s%1s%1s%3s%3s%3s%3s%3s%3s\n",$corx2[$l],$cory2[$l],$corz2[$l],$blanc,$get[3],$blanc,$get[4],$get[5],$get[6],$get[4],$get[5],$get[6];
			}
			else{
				printf OUTC "%10s%10s%10s%1s%2s%3s%3s%3s%3s%3s%3s\n",$corx2[$l],$cory2[$l],$corz2[$l],$blanc,$get[3],$get[4],$get[5],$get[6],$get[4],$get[5],$get[6];
			};
			$rdkatom[$rdki]=$get[3];
                        $rdki++;
		};
	};

## Print Bond

	foreach $l (1..$nbond1){
		@get=split(' ',$lignebond1[$l]);
                               @coller=split(' *',$get[0]);
                                @coller2=split(' *',$get[1]);
                                if(@coller==6 && $get[1] ne ""){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[2]=$get[1];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $get[1] eq ""){
                                        $get[0]=$coller[0].$coller[1];
                                        $get[1]=$coller[2].$coller[3].$coller[4];
                                        $get[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($lignebond1[$l]=~/^\s/){
                                                $get[0]=$coller[0].$coller[1];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($lignebond1[$l]=~/^\s/){
                                                $get[0]=$coller[0];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
						$get[2]=$get[1];
                                                $get[1]=$coller[3];
                                        };					
                                }
                                elsif(@coller2==4){
                                        $get[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $get[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                        $get[2]=$coller[6];
                                };

				$get[3]=0;

		if($get[0] != $iatom1 && $get[1] != $iatom1){
			$get[0]=$get[0]-1 if($iatom1 < $get[0]);
			$get[1]=$get[1]-1 if($iatom1 < $get[1]);
			$get[0]=' '.$get[0] if($get[0] < 100);
			$get[1]=' '.$get[1] if($get[1] < 100);
			printf OUTC "%3s%3s%3s%3s%3s%3s\n",$get[0],$get[1],$get[2],$get[3],$get[3],$get[3]; 
			foreach $w (1..@rdkatom-1){
				if($get[0]==$w){
					$rdkvalence[$w]=$rdkvalence[$w] + $get[2];
				}
				elsif($get[1]==$w){
					$rdkvalence[$w]=$rdkvalence[$w] + $get[2];
				};
			};
		};
	};

	foreach $l (1..$nbond2){ 
		@get=split(' ',$lignebond2[$l]);
                               @coller=split(' *',$get[0]);
                                @coller2=split(' *',$get[1]);
                                if(@coller==6 && $get[1] ne ""){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[2]=$get[1];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $get[1] eq ""){
                                        $get[0]=$coller[0].$coller[1];
                                        $get[1]=$coller[2].$coller[3].$coller[4];
                                        $get[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($lignebond2[$l]=~/^\s/){
                                                $get[0]=$coller[0].$coller[1];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($lignebond2[$l]=~/^\s/){
                                                $get[0]=$coller[0];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
						$get[2]=$get[1];
                                                $get[1]=$coller[3];
                                        };					
                                }
                                elsif(@coller2==4){
                                        $get[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $get[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $get[0]=$coller[0].$coller[1].$coller[2];
                                        $get[1]=$coller[3].$coller[4].$coller[5];
                                        $get[2]=$coller[6];
                                };

				$get[3]=0;

		if($get[0] != $iatom2 && $get[1] != $iatom2){
			$get[0]=$get[0]-1 if($iatom2 < $get[0]);
			$get[0]=$get[0]+$nbatom1-1;
			$get[1]=$get[1]-1 if($iatom2 < $get[1]);
			$get[1]=$get[1]+$nbatom1-1;
			$get[0]=' '.$get[0] if($get[0] < 100);
			$get[1]=' '.$get[1] if($get[1] < 100);
			printf OUTC "%3s%3s%3s%3s%3s%3s\n",$get[0],$get[1],$get[2],$get[3],$get[3],$get[3]; 
			foreach $w (1..@rdkatom-1){
                                if($get[0]==$w){
                                        $rdkvalence[$w]=$rdkvalence[$w] + $get[2];
                                }
                                elsif($get[1]==$w){
                                        $rdkvalence[$w]=$rdkvalence[$w] + $get[2];
                                };
                        };
		};
	}; 	
	$t=1;
	$get[3]=0;
	if($iatom2 > $l2){
		$l2b=$nbatom1-1+$l2;
	}
	else{
		$l2b=$nbatom1-2+$l2;
	};
        if($iatom1 < $l1){
                $l1b=$l1-1;
        }
        else{
                $l1b=$l1;
        };

        $l1b=' '.$l1b if($l1b < 100);
        $l2b=' '.$l2b if($l2b < 100);

	printf OUTC "%3s%3s%3s%3s%3s%3s\n",$l1b,$l2b,$t,$get[3],$get[3],$get[3];
	foreach $w (1..@rdkatom-1){
		if($l1b==$w){
			$rdkvalence[$w]=$rdkvalence[$w] + $t;
                }
		elsif($l2b==$w){
                	$rdkvalence[$w]=$rdkvalence[$w] + $t;
                };
        };
	
	#rdkit 2021
	#M CHG
	$rdkchg="";
	$rdkchgi=0;
	foreach $w (1..@rdkatom-1){
		#print "$w $rdkatom[$w] $rdkvalence[$w]\n";
		if($rdkatom[$w] eq "N" && $rdkvalence[$w]==4){
			print "ici\n";
			$addchg=sprintf " %3s %3s",$w,"1";
			$rdkchg=$rdkchg.$addchg;
			$rdkchgi++;
		};
		if($rdkatom[$w] eq "O" && $rdkvalence[$w]==1){
			$addchg=sprintf " %3s %3s",$w,"-1";
			$rdkchg=$rdkchg.$addchg;
			$rdkchgi++;
		};
	};
	if($rdkchgi>0){
		$rdkchg2=sprintf "M  CHG%3s",$rdkchgi;
		$rdkchg=$rdkchg2.$rdkchg;
		print OUTC "$rdkchg\n";
	};

	#add 17 juin 2008 pb in corina !
	print OUTC "M  END\n";

## Print END

	if($iatom2 > $posr){
		$newposr=$nbatom1-1+$posr;
	}
	else{
		$newposr=$nbatom1-2+$posr;
	};

	#modif 6 mai 2008
	#printf OUTC "> <ID> ($moli)\n";
	#printf OUTC "$moli\n";
	#printf OUTC "\n";
	
	#modif 6 mai 2008
	if($origin != 2){
		printf OUTC "> <FGT>\n";
		printf OUTC "$file1 $noatom1 $file2 $noatom2\n";
		printf OUTC "\n";
	};
		
	$memorylego="";
	$filelego1="";
	if($file1=~/^lego/){
		$filelego1=$file1;
		$filelego1=~s/^lego//;
	};	
	$filelego2="";
	if($file2=~/^lego/){
		$filelego2=$file2;
		$filelego2=~s/^lego//;
	};	
	#print "filelego1=$filelego1 filelego2=$filelego2\n";

		$lignedata3="";

		if($pointsfgt2 ne ''){
			#print "lignedata2 $lignedata2\n"; #only if ORIGIN
	
			@spptfgto=split('-',$lignedata2);
		     	@spptfgt=split('-',$pointsfgt2);
			$spptfgtvu =0;
		     	foreach $sppi (0..@spptfgt-1){
				#print "B 2 $l2 == $spptfgt[$sppi] and $l1 and $spptfgto[$sppi]\n";
				#nov2021
				#if($l2 == $spptfgt[$sppi] && ( $spptfgto[$sppi] eq "$l2_$l1" || $spptfgto[$sppi] eq "$l1_$l2" ) && $spptfgtvu ==0){
				if($l2i==-1 && $l2 == $spptfgt[$sppi] && ( $spptfgto[$sppi] eq "$l2_$l1" || $spptfgto[$sppi] eq "$l1_$l2" ) && $spptfgtvu ==0){
					$spptfgt[$sppi]="";
					$spptfgto[$sppi]="";
					$spptfgtvu =1;
				}
				elsif($sppi==$l2i && $spptfgtvu ==0){
					$spptfgt[$sppi]="";
					$spptfgto[$sppi]="";
					$spptfgtvu =1;
				}
		        	elsif($iatom2 > $spptfgt[$sppi]){
					$memorylego=$memorylego."$filelego2 $spptfgt[$sppi] " if($filelego2 ne "");
					$spptfgt[$sppi]=$nbatom1-1+$spptfgt[$sppi];
					$memorylego=$memorylego."$spptfgt[$sppi]\n" if($filelego2 ne "");
				}
				else{
					$memorylego=$memorylego."$filelego2 $spptfgt[$sppi] " if($filelego2 ne "");
					$spptfgt[$sppi]=$nbatom1-2+$spptfgt[$sppi];
					$memorylego=$memorylego."$spptfgt[$sppi]\n" if($filelego2 ne "");
				};
			};

		        $pointsfgt2=join('-',@spptfgt);
			$pointsfgt2=~s/ //g;
			$pointsfgt2=~s/--/-/g;
			$pointsfgt2=~s/^-//;
			$pointsfgt2=~s/-$//;

			$lignedata2=join('-',@spptfgto);
			$lignedata2=~s/ //g;
			$lignedata2=~s/--/-/g;
			$lignedata2=~s/^-//;
			$lignedata2=~s/-$//;

			if($pointsfgt2 ne ''){
				printf OUTC "> <POINTS>\n";
				printf OUTC "$pointsfgt2";
				#print "$pointsfgt2";
				if($lignedata2 ne ""){
					$lignedata3="> <ORIGIN>\n$lignedata2";
				};	
			};
		};
#print "$lignedata3\n";

		if($pointsfgt1 ne ''){
			#print "lignedata $lignedata\n";
			@spptfgto=split('-',$lignedata);
			@spptfgt=split('-',$pointsfgt1);
			$spptfgtvu =0;
		     	foreach $sppi (0..@spptfgt-1){
				#print "B 1 $l1 == $spptfgt[$sppi]\n";
				#nov2021
				#if($l1 == $spptfgt[$sppi] && ($spptfgto[$sppi] eq "$l1_$l2" || $spptfgto[$sppi] eq "$l2_$l1") && $spptfgtvu ==0){
				if($l1i==-1 && $l1 == $spptfgt[$sppi] && ($spptfgto[$sppi] eq "$l1_$l2" || $spptfgto[$sppi] eq "$l2_$l1") && $spptfgtvu ==0){
					$spptfgt[$sppi]="";
					$spptfgto[$sppi]="";
					$spptfgtvu =1;
				}
				elsif($sppi==$l1i && $spptfgtvu ==0){
					$spptfgt[$sppi]="";
					$spptfgto[$sppi]="";
					$spptfgtvu =1;
				}
		        	elsif($iatom1 < $spptfgt[$sppi]){
					$memorylego=$memorylego."$filelego1 $spptfgt[$sppi] " if($filelego1 ne "");
					$spptfgt[$sppi]=$spptfgt[$sppi]-1;
					$memorylego=$memorylego."$spptfgt[$sppi]\n" if($filelego1 ne "");
				}
				else{
					$memorylego=$memorylego."$filelego1 $spptfgt[$sppi] $spptfgt[$sppi]\n" if($filelego1 ne "");
				};	

			};	

			if($filelego1 eq "" && $origin==2 && $file1 eq "combin.sdf" && $lego1[0] ne ""){
				#print "modif @lego1\n";
				foreach $legoi (0..@lego1-1){
					@splitlego=split(' ',$lego1[$legoi]);
					if($iatom1 < $splitlego[2]){
						$splitlego[2]=$splitlego[2]-1;
						$memorylego=$memorylego."$splitlego[0] $splitlego[1] $splitlego[2]\n";
					}
					else{
						$memorylego=$memorylego."$lego1[$legoi]\n";
					};
				};	
			};
			
			$pointsfgt1=join('-',@spptfgt);
			$pointsfgt1=~s/ //g;
			$pointsfgt1=~s/--/-/g;
			$pointsfgt1=~s/^-//;
			$pointsfgt1=~s/-$//;

                        $lignedata=join('-',@spptfgto);
                        $lignedata=~s/ //g;
                        $lignedata=~s/--/-/g;
                        $lignedata=~s/^-//;
                        $lignedata=~s/-$//;

			if($pointsfgt2 ne '' && $pointsfgt1 ne ''){
				printf OUTC "-$pointsfgt1";
				if($lignedata ne ""){
					$lignedata3=$lignedata3."-$lignedata";
				};	
			}
			elsif($pointsfgt1 ne ''){
				printf OUTC "> <POINTS>\n";
				printf OUTC "$pointsfgt1";
				print "$pointsfgt1";
				$lignedata3="> <ORIGIN>\n$lignedata";
			};
		};

		if($pointsfgt2 ne '' || $pointsfgt1 ne ''){
			printf OUTC "\n\n";
			if($lignedata3 ne ""){
				printf OUTC "$lignedata3\n\n";
			};			
		};

		#modif 6 mai 2008
		if($origin == 2){
			printf OUTC "> <LEGO>\n";
			printf OUTC "$memorylego\n";
		};	
		
		printf OUTC "\$\$\$\$\n";


		print "$file3 DONE\n";

###################################################################
###################################################################
