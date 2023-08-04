#!/usr/bin/perl

        $leaexe=$0;
        #windows requires 2 steps
        $leaexe=~s/properties\.pl$//;
        $leaexe=~s/\/$//;
        #print "perl scripts in $leaexe\n";

	local($filesdf)=@ARGV;	

	if($filesdf eq ""){
		die "usage: .pl <sdf file>\n";	
	};	

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

        %prop=(
                'mm',0.0,
                'logp',0.0,
                'rg',0.0,
                'ix',0.0,
                'iy',0.0,
                'iz',0.0,
                'length',0.0,
                'nbhd',0.0,
                'nbha',0.0,
                'nbatom',0.0,
		'fsp3',0.0,
		'formal_charge',0.0,
        );

	$property=" mm logp rg ix iy iz length nbhd nbha nbatom fsp3 formal_charge";
	#formal_charge can be a zwitterion

	@fprop=split(" ",$property);
	
	#print  "MOLECULAR PROPERTIES :\n";
	open(PPT,">properties.out");
		print PPT "#ID class ";
		foreach $propi (0..@fprop-1){
			print PPT " $fprop[$propi] ";
		};
		print PPT "\n";

		
	##### PROPRIETES FROM SDF
	
	$flagnew=1;
	$moli=0;
	open(MOL,"<$filesdf");
	while(<MOL>){
		if($flagnew){
			$istratom=0;
			$masse=0;
			$flagname=0;
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
			$namegeneric="";
			$namemdl="";
			$zwitt="";
			unlink "tmpi.sdf" if(-e "tmpi.sdf");
			open(OUT,">tmpi.sdf");
                };
                @getstr = split(' ',$_);
		print OUT $_;
		
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

                        #$compt++;
                };
		if($flagname && $getstr[0] ne ''){
			$namemdl=$getstr[0] if($mdl);
			$namegeneric=$getstr[0] if($name);
			$mdl=0;
			$name=0;
			$flagname=0;
		};
		if ($_=~/^>/ && ($_=~/MDLNUMBER/ || $_=~/ZINC/ || $_=~/name/ || $_=~/NAME/)){
			$flagname=1;
			$mdl=1 if($_=~/MDLNUMBER/ || $_=~/ZINC/);
			$name=1 if($_=~/name/ || $_=~/NAME/);
		};	
                if ($_=~/\$\$\$\$/){
			close(OUT);
	
			$f4="tmpi.sdf";
			#chop($tmpmol2 = `$leaexe/SDF_MOL2.pl tmpi.sdf ` );
			#rename "tmpi_1.mol2", "tmpi.mol2";
			#$f3="tmpi.mol2";

			$moli++;	
			
                        $flagnew=1;
			#print "$moli atomes = $istratom\n";

			if($istratom > 0){
				&calprop;
			};	

			if($namemdl ne ""){
                       	 	print PPT " $namemdl $filesdf ";
			}
			else{
				$namegeneric="$moli" if($namegeneric eq "");
				print PPT " $namegeneric $filesdf ";
			};


			#if zwitterion
			if($zwitt ne "" && $zwitt=~/ \+/ && $zwitt=~/ -/){
				$prop{formal_charge}=$setchg.zwitterion;
			}
			else{
				$prop{formal_charge}=$setchg;
			};

                        foreach $propi (0..@fprop-1){
				$prop{$fprop[$propi]}="0.0" if($prop{$fprop[$propi]} eq "");
                        	print PPT " $prop{$fprop[$propi]} ";
                        };
                        print PPT "\n";

			unlink "tmpi.sdf";
			unlink "tmpi.mol2";
			unlink "chargetmpi.mol2";
                };

                if(($_=~/^M CHG/ || $_=~/^M  CHG/) && $setchg==0){
			$zwitt="";
                        @getmoli=split(' ',$_);
                        $countchg=4;
                        foreach $setchgi (4..@getmoli-1){
                                if($setchgi==$countchg && $getmoli[$setchgi] ne ""){
					if($getmoli[$setchgi]=~/-/){
						$zwitt=$zwitt." $getmoli[$setchgi]";
						$getmoli[$setchgi]=~s/-//;
						$setchg=$setchg-$getmoli[$setchgi];
					}
					elsif($getmoli[$setchgi]=~/\+/){
						$zwitt=$zwitt." $getmoli[$setchgi]";
						$getmoli[$setchgi]=~s/\+//;
						$setchg=$setchg+$getmoli[$setchgi];
					}
					else{# no sign +
						$zwitt=$zwitt." +$getmoli[$setchgi]";
						$setchg=$setchg+$getmoli[$setchgi];
                                        };
					$countchg=$countchg+2;
                                };
                        };
			#$setchg=$setchg."zwitterion" if($zwitt=~/ \+/ && $zwitt=~/ -/ && $setchg==0);
                };
        };
        close(MOL);

        close(PPT);



################################################################################
#
sub calprop{

	$exemm=1;
	$exelogp=1;
	$exerg=1;
	$exeix=1;
	$exeiy=1;
	$exeiz=1;
	$exelength=1;
	$exenbhd=1;
	$exenbha=1;
	$exenbatom=1;
	$exefsp3=1;
	
#================================= Masse
#       $masse
if($exemm || $exelipinski){
        $prop{mm}=sprintf"%4.2f",$masse;
};

#================================= Nbatom
if($exenbatom || $exelipinski){
        $prop{nbatom}=$atomlourd;
};

#================================= Rayon de giration

if($exeix || $exeiy || $exeiz || $exerg){
        $rgx=0;
        $rgy=0;
        $rgz=0;
        foreach $ig (1..$istratom){
                $rgx=$rgx+$strx[$ig]*$tabmm{$atom[$ig]};
                $rgy=$rgy+$stry[$ig]*$tabmm{$atom[$ig]};
                $rgz=$rgz+$strz[$ig]*$tabmm{$atom[$ig]};
        };

#center of mass:
        $rgx=$rgx/$masse;
        $rgy=$rgy/$masse;
        $rgz=$rgz/$masse;

        $rayongyration = 0;
        $raygx=0;
        $raygy=0;
        $raygz=0;
        $rayg=0;

        foreach $ig (1..$istratom){
                $raygx=$rgx-$strx[$ig];
                $raygx=$raygx*$raygx;
                $raygy=$rgy-$stry[$ig];
                $raygy=$raygy*$raygy;
                $raygz=$rgz-$strz[$ig];
                $raygz=$raygz*$raygz;
                $rayg=sqrt($raygx+$raygy+$raygz);
                $rayg=$rayg*$rayg;
                $rayongyration = $rayongyration + $rayg;
        };
        $rayongyration = (sqrt($rayongyration))/($istratom+1);
        $prop{rg}=sprintf"%4.2f",$rayongyration;
};

#================================= Length and Axes Inertie Ix, Iy, Iz

$oldfashion=0;

if($exelength && $oldfashion==0){
#length
	$biglength=-1;
	$strlengthimin="";
	$strlengthimax="";
	foreach $ig (1..$istratom){
		foreach $im ($ig+1..$istratom){
			$length=sqrt(($strx[$ig]-$strx[$im])**2 + ($stry[$ig]-$stry[$im])**2 + ($strz[$ig]-$strz[$im])**2);
			if($length > $biglength){
				$biglength=$length;
				$strlengthimin=$ig;
				$strlengthimax=$im;
			};
		};	
	};
	#print "length $strlengthimin -> $strlengthimax = $biglength\n";
	$prop{length}=sprintf"%4.2f",$biglength;
};


if(($exeix || $exeiy || $exeiz) && $oldfashion==0){

#diagonalized components Ixx, Iyy and Izz of the moment of inertia tensor

	#center of mass is calculated above ($rgx, $rgy, $rgz)
	$ixx=0;
	$iyy=0;
	$izz=0;
	foreach $ig (1..$istratom){
		$ixx=$tabmm{$atom[$ig]}*(($rgy-$stry[$ig])*($rgy-$stry[$ig]) + ($rgz-$strz[$ig])*($rgz-$strz[$ig]));
		$iyy=$tabmm{$atom[$ig]}*(($rgx-$strx[$ig])*($rgx-$strx[$ig]) + ($rgz-$strz[$ig])*($rgz-$strz[$ig]));
		$izz=$tabmm{$atom[$ig]}*(($rgx-$strx[$ig])*($rgx-$strx[$ig]) + ($rgy-$stry[$ig])*($rgy-$stry[$ig]));
	};

	if($ixx < $iyy){
		$itmp=$iyy;
		$iyy=$ixx;
		$ixx=$itmp;
	};
	if($ixx < $izz){
		$itmp=$izz;
		$izz=$ixx;
		$ixx=$itmp;
	};
	if($iyy < $izz){
		$itmp=$iyy;
		$iyy=$izz;
		$izz=$itmp;
	};

	#inverse for mimicking paper of akritopoulou-zanze et al, DDT, 12, 948-952
	$itmp=$ixx;	
	$ixx=$izz;
	$izz=$itmp;

	#Normalized for mimicking paper of akritopoulou-zanze et al, DDT, 12, 948-952
	$prop{ix}=sprintf"%4.2f",$ixx/$izz if($izz!=0);
	$prop{iy}=sprintf"%4.2f",$iyy/$izz if($izz!=0);
	$prop{iz}=sprintf"%4.2f",$izz/$izz if($izz!=0);

}
elsif(($exeix || $exeiy || $exeiz || $exelength) && $oldfashion==1){

        $strlongxmin=$strx[1];
        $strlongxmax=$strx[1];
        $strlongymin=$stry[1];
        $strlongymax=$stry[1];
        $strlongzmin=$strz[1];
        $strlongzmax=$strz[1];

        $strlengthiminx=1;
        $strlengthimaxx=1;
        $strlengthiminy=1;
        $strlengthimaxy=1;
        $strlengthiminz=1;
        $strlengthimaxz=1;

        foreach $longxk (2..$istratom){
                if ($atom[$longxk] ne 'H'){
                        if ($strx[$longxk] > $strlongxmax){
                                $strlongxmax=$strx[$longxk];
                                $strlengthimaxx=$longxk;
                        };
                        if ($strx[$longxk] < $strlongxmin){
                                $strlongxmin=$strx[$longxk];
                                $strlengthiminx=$longxk;
                        };
                        if ($stry[$longxk] > $strlongymax){
                                $strlongymax=$stry[$longxk];
                                $strlengthimaxy=$longxk;
                        };
                        if ($stry[$longxk] < $strlongymin){
                                $strlongymin=$stry[$longxk];
                                $strlengthiminy=$longxk;
                        };
                        if ($strz[$longxk] > $strlongzmax){
                                $strlongzmax=$strz[$longxk];
                                $strlengthimaxz=$longxk;
                        };
                        if ($strz[$longxk] < $strlongzmin){
                                $strlongzmin=$strz[$longxk];
                                $strlengthiminz=$longxk;
                        };
                };
        };
        $longueurx=$strlongxmax - $strlongxmin;
        $longueurx=$longueurx*$longueurx;
        $longueurx=sqrt($longueurx);

        $longueury=$strlongymax - $strlongymin;
        $longueury=$longueury*$longueury;
        $longueury=sqrt($longueury);

        $longueurz=$strlongzmax - $strlongzmin;
        $longueurz=$longueurz*$longueurz;
        $longueurz=sqrt($longueurz);

        if (($longueurx >= $longueury) && ($longueury >= $longueurz)){
                $ix=$longueurx;
                $iy=$longueury;
                $iz=$longueurz;
                $strlengthimin=$strlengthiminx;
                $strlengthimax=$strlengthimaxx;
        }
        elsif (($longueurx >= $longueurz) && ($longueurz >= $longueury)){
                $ix=$longueurx;
                $iy=$longueurz;
                $iz=$longueury;
                $strlengthimin=$strlengthiminx;
                $strlengthimax=$strlengthimaxx;
        }
        elsif (($longueury >= $longueurx) && ($longueurx >= $longueurz)){
                $ix=$longueury;
                $iy=$longueurx;
                $iz=$longueurz;
                $strlengthimin=$strlengthiminy;
                $strlengthimax=$strlengthimaxy;
        }
        elsif (($longueury >= $longueurz) && ($longueurz >= $longueurx)){
                $ix=$longueury;
                $iy=$longueurz;
                $iz=$longueurx;
                $strlengthimin=$strlengthiminy;
                $strlengthimax=$strlengthimaxy;
        }
        elsif (($longueurz >= $longueurx) && ($longueurx >= $longueury)){
                $ix=$longueurz;
                $iy=$longueurx;
                $iz=$longueury;
                $strlengthimin=$strlengthiminz;
                $strlengthimax=$strlengthimaxz;
        }
       elsif (($longueurz >= $longueury) && ($longueury >= $longueurx)){
                $ix=$longueurz;
                $iy=$longueury;
                $iz=$longueurx;
                $strlengthimin=$strlengthiminz;
                $strlengthimax=$strlengthimaxz;
        };

	#length is defined by the atoms used to calculate the biggest inertial axis (and not the 2 atoms that are the most far away from each other)
	
        $lengthix=($strx[$strlengthimax]-$strx[$strlengthimin]);
        $lengthiy=($stry[$strlengthimax]-$stry[$strlengthimin]);
        $lengthiz=($strz[$strlengthimax]-$strz[$strlengthimin]);
        $lengthix=$lengthix*$lengthix;
        $lengthiy=$lengthiy*$lengthiy;
        $lengthiz=$lengthiz*$lengthiz;
        $length = $lengthix + $lengthiy + $lengthiz;
        $length = sqrt($length);

        $prop{length}=sprintf"%4.2f",$length;
        $prop{ix}=sprintf"%4.2f",$ix;
        $prop{iy}=sprintf"%4.2f",$iy;
        $prop{iz}=sprintf"%4.2f",$iz;
};

#================================= Nombre de Ha et Hd
if($exenbha || $exenbhd || $exelipinski){

        $nbha=0;
        $nbhd=0;
        foreach $ig (1..$istratom){
                $nbhd++ if ($atom[$ig] eq 'H' && ($fonc[$ig]=~/ 1-O /||$fonc[$ig]=~/ 1-N /));
                $nbha++ if ($atom[$ig] eq 'N' && $coval[$ig] <= 3);
                $nbha++ if ($atom[$ig] eq 'O');
        };

        $prop{nbha}=$nbha;
        $prop{nbhd}=$nbhd;
};

#================================= Fsp3 = nbCsp3 / nbC
if($exefsp3){
	$fsp3=0;
	$nbc=0;
	foreach $ig (1..$istratom){
		$nbc++ if ($atom[$ig] eq 'C');
		@getfs=split(' ',$covfonc[$ig]);
		$nbvoisin=@getfs;
		$fsp3++ if($nbvoisin==4 && $atom[$ig] eq 'C');
	};	
	if($nbc>0){
	$prop{fsp3}=sprintf"%4.2f",($fsp3/$nbc);
	}
	else{
	$prop{fsp3}=0;
	};
	#print "Csp3 = $fsp3 ; Ctot = $nbc ; fsp3 = $prop{fsp3}\n";
};

#================================= LogP by XLOGP (atomic contribution by Wang et al.)
if($exelogp || $exelipinski){

#now uses rdkit MolLogP	

	if(-e $f4){
		chop($tmplog=`python $leaexe/rdkit-MolLogP.py $f4`);
		#chop($tmplog=` $leaexe/rdkit-MolLogP.py $f4`);
		@getlogp=split(' ',$tmplog);
		if($getlogp[1] ne ""){
			$prop{logp}=$getlogp[1];
		}
		else{
			$prop{logp}=0.0;
		};

        }
        else{
                $prop{logp}=0.0;
        };
};

};
