#!/usr/bin/perl

# UPGRADES that will change the KEY files of databases
# 18 april 2011: detect phosphines R/S
# 18 april 2011: detect imines -N=C()()  trans/cis
# Don't take account for symetry like mesoR or mesoS  
# accept bore atoms December 2016
# -CO2 corrected Jan2017	
# Pb in recursivity in C*: refeuille3 corrected

# FOR DEBUGGING : set debug=1 and debog=1

#Keep the following print as the reply to the require perl function
print "";

sub keymol2{

	local($fileopen,$filekey,$eltx,$flagframework)=@_;

        $leaexe=$0;
        #windows requires 2 steps
	#called by main program CLASS_FGT.pl
        $leaexe=~s/lea3d-CLASS_FGT\.pl$//;
        $leaexe=~s/\/$//;
        #print "perl scripts in $leaexe\n";
	$cyclesdf=$leaexe."/lea3d-CYCLE_SDF.pl";
	require $cyclesdf;

	# create a mode which allows the framework evaluation
	# single bond, carbons and no hydrogens
	$framework=0;
	$framework=1 if($flagframework);
	$halo=" F Cl Br I H X ";

	$cluster=0;
	#$cluster=1 if($flagframework && -e "$filekey" && !-z "$filekey"); # create clusters of framework if $filekey exists 
	
	$fileopen2=$fileopen;
	
	@key_a='';
	@key_b='';
	@occurence='';
	@origin='';
	$flagnewfilekey=-1;
	if(-e "$filekey" && !-z "$filekey"){
		$i=0;
		open(IN,"<$filekey");
		while(<IN>){
		
			@get2=split(' ',$_);
			$origin[$i]="$get2[1] $get2[2]";
		
			@get=split('\*',$get2[0]);
			$key_a[$i]=$get[0];
			$get[0]='';
			$key_b[$i]=join('*',@get);
			$key_b[$i]=~s/^\*//;
			$key_b[$i]=~s/\n$//;
			#$occurence[$i]=1;
			$occurence[$i]=0;
			#print "$key_a[$i] \n$key_b[$i]\n";
			$i++;
			$flagnewfilekey++;
		};
		close(IN);
		$filekey=~s/(.*)\/(.*)/$2/;
		$filekey="additional_".$filekey;
		#print "$filekey created\n";
	}
	else{
		#print "new $filekey\n";
		#system("touch $filekey");
		open(OUT,">$filekey");
		close(OUT);
	};
							
	
	unlink "same.sdf" if(-e "same.sdf");
	unlink "same_additional_key.sdf" if(-e "same_additional_key.sdf");
	unlink "diff.sdf" if(-e "diff.sdf");
	unlink "exclude.sdf" if(-e "exclude.sdf");
	
	$nombre_molecule=0;
	$flagnew=1;
	open(MOL,"<$fileopen");
	while(<MOL>){
		if($flagnew){
			$masse=0;
			$getinfo=0;
			$info='';
			$compt=0;
			$ig=1;
			$jg=0;
			@strx='';
			@stry='';
			@strz='';
			@atom='';
			@atomx='';
			@coval='';
			@fonc='';
			@ifonc='';
			@covfonc='';
			@bond='';
			@listb='';
			@typeb='';
			$blanc=' ';	
			@radius='';
			@lignebond='';	
			@ligne='';	
			$atomlourd=0;
			$flagnew=0;
			$nombre_molecule++;
			$same=0;
			$weird=0;
		};
		@getstr = split(' ',$_);
		$ligne[$compt]=$_;
		$compt++;
		
		if($getinfo){
			$getinfo=0;
			$jointxt=join(' ',@getstr);
			#$info=$info.":".$jointxt;	
			$info=$getstr[0];
		};
		
		if (($compt > 4) && ($ig <= $istratom)){
			$strx[$ig]=$getstr[0];
			$stry[$ig]=$getstr[1];
			$strz[$ig]=$getstr[2];
			
			$atom[$ig]=$getstr[3];
			if($getstr[3] eq "X"){
				$getstr[3]="H";
				$atom[$ig]="H";
				$atomx[$ig]="X";
			}
			else{
				 $atomx[$ig]="-";
			};
			
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
				
				$bond[$getstr[0]]=$bond[$getstr[0]].$blanc.$getstr[1].$blanc.$getstr[2];
				$listb[$getstr[0]]=$listb[$getstr[0]].$blanc.$getstr[1];
				$typeb[$getstr[0]]=$typeb[$getstr[0]].$blanc.$getstr[2];

				$bond[$getstr[1]]=$bond[$getstr[1]].$blanc.$getstr[0].$blanc.$getstr[2];
				$listb[$getstr[1]]=$listb[$getstr[1]].$blanc.$getstr[0];
				$typeb[$getstr[1]]=$typeb[$getstr[1]].$blanc.$getstr[2];


				$fonc[$getstr[0]]=$fonc[$getstr[0]].$blanc.$getstr[2].'-'.$atom[$getstr[1]].$blanc;
				$ifonc[$getstr[0]]=$ifonc[$getstr[0]].$blanc.$getstr[1].$blanc;
				$covfonc[$getstr[0]]=$covfonc[$getstr[0]].$blanc.$getstr[2];
				$coval[$getstr[0]]=$coval[$getstr[0]]+$getstr[2];

				$fonc[$getstr[1]]=$fonc[$getstr[1]].$blanc.$getstr[2].'-'.$atom[$getstr[0]].$blanc;
				$ifonc[$getstr[1]]=$ifonc[$getstr[1]].$blanc.$getstr[0].$blanc;
				$covfonc[$getstr[1]]=$covfonc[$getstr[1]].$blanc.$getstr[2];
				$coval[$getstr[1]]=$coval[$getstr[1]]+$getstr[2];
				$lignebond[$jg]=$_;
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
		#print "$istratom $istrbond\n";

			$flagnew=1;
			$prop{atomlourd}=$atomlourd;
			&convert;
			#print "Mol$moli $mol[$moli-1] converted\n" if($param{VERBOSITY} >= 1);
			
			if($weird){
				open(OUT,">>exclude.sdf");
				foreach $i (0..@ligne-1){
					print OUT $ligne[$i];
				};
				close(OUT);
			}
			elsif($same){
				if($samei > $flagnewfilekey && $flagnewfilekey > -1){
					open(OUT,">>same_additional_key.sdf");
					foreach $i (0..@ligne-1){
						print OUT $ligne[$i];
					};
					close(OUT);
				}
				else{
					open(OUT,">>same.sdf");	
					foreach $i (0..@ligne-1){
						print OUT $ligne[$i];
					};
					close(OUT);
					if($cluster){
						$bi2=$samei+1;
						open(OUT,">>cluster_$bi2.sdf");
						foreach $i (0..@ligne-1){
							print OUT $ligne[$i];
						};
						close(OUT);
					};	
				};	
			}
			else{	
				open(OUT,">>diff.sdf");
				foreach $i (0..@ligne-1){
					print OUT $ligne[$i];
				};
				close(OUT);
			};
		};
		
		if ($_=~/>/ && $_=~/<ID>/){
			@getmoli=split(' ',$_);
			$moli=$getmoli[2];
			$moli=~s/\(//;
			$moli=~s/\)//;
		};
		
		#if ($_=~/>/ && ($_=~/CAS/ || $_=~/ID/ || $_=~/ENTRY/ || $_=~/code_number/ || $_=~/REGNO/ || $_=~/COMP_NAME/ || $_=~/NOM/ || $_=~/NAME/)){
		#	    $getinfo=1;
		#    $info=$info." $getstr[1] $getstr[2]";
		#};
		
		if ($_=~/>/ && $_=~/MDLNUMBER/){
			$getinfo=1;	
		};	
	};
	close(MOL);

	print ('*'x80);
	print "\nkey_number\tOriginal_molecule\t\tOccurence_in_current_file\n"; 
	foreach $bi (0..@key_a-1){
		$bi2=$bi+1;
		if($bi > $flagnewfilekey){
			print "\tkeynumber_add $bi2\t$origin[$bi]\t\t$occurence[$bi]\n";
		}
		else{
			print "\tkeynumber $bi2\t$origin[$bi]\t\t$occurence[$bi]\n";
		};	
	};
	print ('*'x80);
	print "\n";
	
};


###################################################################################
###################################################################################


sub convert{

	$exclude_framework="";
	@type='';
	foreach $bi (1..$istratom){
#		print"$bi _ $atom[$bi] $fonc[$bi] \n";
#		print"$bi _ $atom[$bi] $ifonc[$bi] \n";

		if($atom[$bi] eq 'C'){
			$car=$fonc[$bi];
			$nocc=$car=~s/ 2-/ 2-/g;
			if($fonc[$bi] =~/ 2-/){
				$type[$bi]="C.2";
			}
			elsif($fonc[$bi] =~/ 3-/){
				$type[$bi]="C.1";
			}
			elsif($nocc==2){
				$type[$bi]="C.1";
			}
			elsif($fonc[$bi] =~/ 4-/){
				$type[$bi]="C.ar";	
			}
			else{
				$type[$bi]="C.3";
			};
		};
		
		if($atom[$bi] eq 'N'){
			if($fonc[$bi] =~/ 2-/){
				$type[$bi]="N.2";
			}
			elsif($fonc[$bi] =~/ 3-/){
				$type[$bi]="N.1";
			}
			elsif($fonc[$bi] =~/ 4-/){
				$type[$bi]="N.ar";	
			}
			else{
				$type[$bi]="N.3";
			};
			$type[$bi]="N.4" if($coval[$bi] == 4 && $fonc[$bi] !~/ 2-/ && $fonc[$bi] !~/ 3-/);
			
			$type[$bi]="N.pl3" if($coval[$bi] == 5 && @fonc[$bi]==3); #si NO2
		};
		
		if($atom[$bi] eq 'O'){
		#print "fonc $fonc[$bi]\n";
			if($fonc[$bi] =~/ 2-/){
				$type[$bi]="O.2";
			}
			else{
				$type[$bi]="O.3";
			};
		};
		if($atom[$bi] eq 'S'){
		#print "fonc $fonc[$bi]\n";
			if($fonc[$bi] =~/ 2-/){
				$type[$bi]="S.2";
			}
			else{
				$type[$bi]="S.3";
			};
			if($fonc[$bi] =~/ 2-O /){
				$car=$fonc[$bi];
				$nocc=$car=~s/ 2-O / 2-O /g;
				$type[$bi]="S.o" if($nocc==1);
				$type[$bi]="S.o2" if($nocc==2);
			}
		};
		if($atom[$bi] eq 'P'){
			$type[$bi]="P.3";
		};
		if($atom[$bi] eq 'H'){
			$type[$bi]="H";
		};
		if($atom[$bi] eq 'F'){
			$type[$bi]="F";
		};
		if($atom[$bi] eq 'Cl'){
			$type[$bi]="Cl";
		};
		if($atom[$bi] eq 'Br'){
			$type[$bi]="Br";
		};
		if($atom[$bi] eq 'I'){
			$type[$bi]="I";
		};
		if($atom[$bi] eq 'Si'){
			$type[$bi]="Si";
		};
#Add 14 december 2016
		if($atom[$bi] eq 'B'){
			$type[$bi]="B";
		};

		#special cases
		if($type[$bi] eq ""){
			print "special atom $bi $atom[$bi] not recognized !\n";
			$type[$bi]=$atom[$bi];
		};
	};

	foreach $bi (1..$istratom){
		if($atom[$bi] eq 'C'){

# Fonction acide (O.co2)
			if($fonc[$bi] =~/ 2-O / && $fonc[$bi] =~/ 1-O / ){
			$car=$fonc[$bi];
			$nocc=$car=~s/ 1-O / 1-O /g;
			if ($nocc == 1){ 
				@det=split(' ',$fonc[$bi]);
				$p=0;
				foreach $k (0..@det-1){
					$p=$k if($det[$k] =~/1-O/);
				};
				@det1=split(' ',$ifonc[$bi]);
				$iatom=$det1[$p];
				@det2=split(' ',$ifonc[$iatom]);

# if ($fonc[$iatom] =~/ 1-H /) O.3 et O.2 sinon O.co2 sur les 2 O
				if(@det2==1){
					$type[$iatom]="O.co2";
					$p=0;
					foreach $k (0..@det-1){
						$p=$k if($det[$k] =~/2-O/);
					};
					@det1=split(' ',$ifonc[$bi]);
					$iatom1=$det1[$p];
					$type[$iatom1]="O.co2";
				};
			};
			};
			
# C.cat in guanidinium
			if($fonc[$bi] =~/ 2-N / && $fonc[$bi] =~/ 1-N / ){
				$car=$fonc[$bi];
				$nocc=$car=~s/ 1-N / 1-N /g;
				if($nocc == 2){
					@det=split(' ',$fonc[$bi]);
					@det1=split(' ',$ifonc[$bi]);
					foreach $k (0..@det-1){
						if($det[$k] =~/2-N/){
							$iatom=$det1[$k];
							$car2=$fonc[$iatom];
							$nocc2=$car2=~s/ 1-H / 1-H /g;
							if($nocc2==2){ #means N is protonated
								$type[$bi]="C.cat";
							};
						};
					};
				};
			};	

# Fonction amide (N.am)
			if($fonc[$bi] =~/ 2-O / && $fonc[$bi] =~/ 1-N / ){
				@det=split(' ',$fonc[$bi]);
				foreach $k (0..@det-1){
					if($det[$k] =~/1-N/){
						@det1=split(' ',$ifonc[$bi]);
						$type[$det1[$k]]="N.am" if($atom[$det1[$k]] eq 'N');
					};
				};
			};
		};
		
		if($atom[$bi] eq 'P'){
		
# Fonction phosphate une seule charge - (P et O.co2)

			if($fonc[$bi] =~/ 2-O / && $fonc[$bi] =~/ 1-O / ){
			$car=$fonc[$bi];
			$nocc=$car=~s/ 1-O / 1-O /g;
			$dejavu=0;
			$p=-1;
			if ($nocc == 3){
				@det=split(' ',$fonc[$bi]);
				@det1=split(' ',$ifonc[$bi]);
				foreach $k (0..@det-1){
					if($det[$k] =~/1-O/){
						@det2=split(' ',$ifonc[$det1[$k]]);
						if(@det2==1){
							$type[$det1[$k]]="O.co2";
							$p=$det1[$k];
						};
						$dejavu++ if($fonc[$det1[$k]] =~/ 1-H /);
					};
				};
				
				#if($dejavu>=1 && $p != -1){
				if($p != -1){                     # ici met o.co2 sur tout atom ionis�et =O du phosphate
				        foreach $k (0..@det-1){
					if($det[$k] =~/2-O/){
						@det2=split(' ',$ifonc[$det1[$k]]);
						if(@det2==1){
							$type[$det1[$k]]="O.co2";
							$type[$p]="O.co2";
						};
					};
					};
				};
				
			};
			};			
              };		

	};


				
	
# AROMATIQUES N.ar et C.ar

		#print"Mol$moli\n";
		$nbatom=$istratom;
		&cyclesdf2;
		
if(($cyclef[0] ne "" && $framework) || $framework==0){# if acyclic then skipped conversion !

# N.pl3

		foreach $bi (1..$istratom){
			if($atom[$bi] eq 'N' && $type[$bi] eq 'N.3'){
				
				@det=split(' ',$ifonc[$bi]);
				$p=0;
				foreach $k (0..@det-1){
					$p++ if($type[$det[$k]] eq 'C.cat');
					$p++ if($type[$det[$k]] eq 'C.ar' || $type[$det[$k]] eq 'C.2');
					$p++ if($type[$det[$k]] eq 'C.1');
					$p++ if($type[$det[$k]] eq 'N.2' || $type[$det[$k]] eq 'N.ar'|| $type[$det[$k]] eq 'S.o' || $type[$det[$k]] eq 'S.o2');
				};
				
				$car=$fonc[$bi];
				$nocc=$car=~s/ 1-H / 1-H /g;
				$noccb=$car=~s/ 1-C / 1-C /g;
				
				$noncycle=1;
				$sicycle=0;
				foreach $lc (0..@cyclef-1){
					$tp=' '.$cyclef[$lc].' ';
					$sicycle=1  if($tp=~/ $bi /);
					$noncycle=0  if(($tp=~/ $bi /) && $cyclefar[$lc] ne "ar");
				};
				#print "$p $sicycle $noncycle $noccb\n";
				$type[$bi]="N.pl3" if($p >=1 && $sicycle==1 && $noncycle==0 && $nocc==1 && $type[$bi] eq "N.ar");

				$type[$bi]="N.pl3" if($p >= 2);
				$type[$bi]="N.pl3" if($p==1 && $nocc==2 && $sicycle==0);
				$type[$bi]="N.pl3" if($p==1 && $noccb>=1 && $nocc==1 && $sicycle==0);
				$type[$bi]="N.pl3" if($p==1 && $noccb>=2 && $sicycle==0);
				
				$type[$bi]="N.pl3" if($p >=1 && $sicycle==1 && $noncycle==1);
				
				#$type[$bi]="N.pl3" if($p >=1 && $sicycle==1 && $noncycle==0 && $nocc==1 && @arom==6);
				#$type[$bi]="N.pl3" if($p >=1 && $sicycle==1 && $noncycle==0 && $nocc==1 && @arom==5);
				#$type[$bi]="N.pl3" if($p >=1 && $sicycle==1 && $noncycle==0 && $noccb>=2);
				
				
			}  #si NO2
			elsif($atom[$bi] eq 'N' && $type[$bi] eq 'N.2'){
				$car=$fonc[$bi];
			        $nocc=$car=~s/ 2-O / 2-O /g;
			        $noccb=$car=~s/ 1-O / 1-O /g;
			        $type[$bi]="N.pl3" if($nocc == 1 && $noccb== 1);
			
			        @det=split(' ',$fonc[$bi]);
	                        foreach $k (0..@det-1){
	                        	if($det[$k] =~/1-O/){
	                        		@det1=split(' ',$ifonc[$bi]);
	                        		$type[$det1[$k]]="O.2" if($atom[$det1[$k]] eq 'O');
	                        	};
	                        };
			};
		};

# S.2 si S.3 pres de C.ar ou C.2

		foreach $bi (1..$istratom){
			if($atom[$bi] eq 'S' && $type[$bi] eq 'S.3'){
				@det=split(' ',$ifonc[$bi]);
				foreach $k (0..@det-1){
					$type[$bi]="S.2" if($type[$det[$k]] eq 'C.ar' || $type[$det[$k]] eq 'C.2' || $type[$det[$k]] eq 'C.1' || $type[$det[$k]] eq 'N.ar' || $type[$det[$k]] eq 'N.2' );
				};
			};
		};


### CONVERSION
		
	#Variables en vue du classement des FGTS
	
	$class_c3=0;
	$class_c2=0;
	$class_c1=0;
	$class_car=0;
	$class_cat=0;
	$class_n3=0;
	$class_n2=0;
	$class_n1=0;
        $class_nar=0;
        $class_nam=0;
        $class_npl3=0;
        $class_n4=0;
        $class_o3=0;
        $class_o2=0;
        $class_oco2=0;

        $class_ospc=0;
        $class_ot3p=0;

        $class_s3=0;
        $class_s2=0;
        $class_so=0;
        $class_so2=0;
        $class_p3=0;
        $class_h=0;

        $class_hspc=0;
        $class_ht3p=0;

        $class_f=0;
        $class_cl=0;
        $class_br=0;
        $class_i=0;
        $class_si=0;

        $class_lp=0;
        $class_du=0;
        $class_na=0;
        $class_k=0;
        $class_ca=0;
        $class_li=0;
        $class_al=0;


	$blanc=' ';
	$f4="mol$moli.mol2";
	#open(OUTC,">$f4");
	
		if($info ne ''){
			#print OUTC "#	  Creating by LEA \n";
			#print OUTC "#  $info \n";
			#print OUTC "#\n";
		}
		else{
			#print OUTC "#\n";
			#print OUTC "#	  Creating by LEA \n";
			#print OUTC "#\n";
		};
		
		#print OUTC "\n";
		#print OUTC "@<TRIPOS>MOLECULE\n";
		#print OUTC "mol$moli\n";
		#printf OUTC "%4s%1s%4s\n",$istratom,$blanc,$istrbond;
		#print OUTC "SMALL\n";
		#print OUTC "NO_CHARGES\n\n\n";
		#print OUTC "@<TRIPOS>ATOM\n";
		foreach $bi (1..$istratom){
			$a1=$bi;
			$a2="$atom[$bi]$bi";
			$a3=sprintf "%7.4f",$strx[$bi];
			$a4=sprintf "%7.4f",$stry[$bi];
			$a5=sprintf "%7.4f",$strz[$bi];
			$a6=$type[$bi];
			#printf OUTC "%4s%5s%6s%11s%11s%11s $a6\n",$a1,$a2,$blanc,$a3,$a4,$a5;
			
			if($framework){
				$typebiold=$type[$bi];
				if($halo=~/ $atom[$bi] /){
					$type[$bi]="H";
				}
				else{	
#print "$bi \n";
					$type[$bi]="C.3";

					# check if terminal C=O or terminal C=C(X)X or terminal C=N(X)
					$sicycle=0;
                                	foreach $lc (0..@cyclef-1){
                                        	$tp=' '.$cyclef[$lc].' ';
                                        	$sicycle=1  if($tp=~/ $bi /);
                                	};
					if($sicycle==0){
						@det=split(' ',$fonc[$bi]);
						if(@det==1 && $fonc[$bi]=~/2-/ && ($atom[$bi] eq 'S' || $atom[$bi] eq 'O')){
							$exclude_framework=$exclude_framework." $bi ";
							$type[$bi]="H";
						}
						elsif(@det==2 && $atom[$bi] eq 'N' && $fonc[$bi]=~/2-/){
							@det2=split(' ',$ifonc[$bi]);
							$nbfx=0;
							foreach $k (0..@det-1){
								$nbfx++ if($det[$k]=~/1-/ && $atomx[$det2[$k]] eq 'X');
							};
							if($nbfx==1){
								$exclude_framework=$exclude_framework." $bi ";
								$type[$bi]="H";
							};
						}
						elsif(@det==3 && $atom[$bi] eq 'C'){
							@det2=split(' ',$ifonc[$bi]);
							$nbfx=0;
							foreach $k (0..@det-1){
								$nbfx++ if($det[$k]=~/1-/ && $atomx[$det2[$k]] eq 'X');
							};
							if($nbfx==2){
								$exclude_framework=$exclude_framework." $bi ";
								$type[$bi]="H";
							};
						};
					};
				};
			};
#print "$exclude_framework\n";
	#	print "$bi \n" if($type[$bi]=~/C\.3/);	
			$class_c3++ if($type[$bi]=~/C\.3/);
			$class_c2++ if($type[$bi]=~/C\.2/);
			$class_c1++ if($type[$bi]=~/C\.1/);
			$class_car++ if($type[$bi]=~/C\.ar/);
			$class_cat++ if($type[$bi]=~/C\.cat/);
			$class_n3++ if($type[$bi]=~/N\.3/);
			$class_n2++ if($type[$bi]=~/N\.2/);
			$class_n1++ if($type[$bi]=~/N\.1/);
       	 		$class_nar++ if($type[$bi]=~/N\.ar/);
        		$class_nam++ if($type[$bi]=~/N\.am/);
        		$class_npl3++ if($type[$bi]=~/N\.pl3/);
        		$class_n4++ if($type[$bi]=~/N\.4/);
        		$class_o3++ if($type[$bi]=~/O\.3/);
        		$class_o2++ if($type[$bi]=~/O\.2/);
        		$class_oco2++ if($type[$bi]=~/O\.co2/);

       		 	$class_ospc++ if($type[$bi]=~/O\.spc/);
        		$class_ot3p++ if($type[$bi]=~/O\.t3p/);

        		$class_s3++ if($type[$bi]=~/S\.3/);
        		$class_s2++ if($type[$bi]=~/S\.2/);
        		$class_so++ if($type[$bi]=~/S\.o/);
        		$class_so2++ if($type[$bi]=~/S\.o2/);
        		$class_p3++ if($type[$bi]=~/P\.3/);
        		$class_h++ if($type[$bi]=~/H/);

        		$class_hspc++ if($type[$bi]=~/H\.spc/);
        		$class_ht3p++ if($type[$bi]=~/H\.t3p/);

       	 		$class_f++ if($type[$bi]=~/F/);
        		$class_cl++ if($type[$bi]=~/Cl/);
        		$class_br++ if($type[$bi]=~/Br/);
        		$class_i++ if($type[$bi]=~/I/);
        		$class_si++ if($type[$bi]=~/Si/);

        		$class_lp++ if($type[$bi]=~/LP/);
        		$class_du++ if($type[$bi]=~/Du/);
        		$class_na++ if($type[$bi]=~/Na/);
       	 		$class_k++ if($type[$bi]=~/K/);
        		$class_ca++ if($type[$bi]=~/Ca/);
        		$class_li++ if($type[$bi]=~/Li/);
        		$class_al++ if($type[$bi]=~/Al/);
		
			if($framework){
				$type[$bi]=$typebiold;
				$class_h=0;
			};
			
		};
	
		$class_ar=0;
		$class_1=0;
		$class_2=0;
		$class_3=0;
		$class_am=0;
	
		@typebondmol2='';
		#print OUTC "@<TRIPOS>BOND\n";
		foreach $bi (1..$istrbond){

			@extract=split(' ',$lignebond[$bi]);

                               @coller=split(' *',$extract[0]);
                                @coller2=split(' *',$extract[1]);
                                if(@coller==6 && $extract[1] ne ""){
                                        $extract[0]=$coller[0].$coller[1].$coller[2];
                                        $extract[2]=$extract[1];
                                        $extract[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $extract[1] eq ""){
                                        $extract[0]=$coller[0].$coller[1];
                                        $extract[1]=$coller[2].$coller[3].$coller[4];
                                        $extract[2]=$extract[5];
                                }
                                elsif(@coller==5){
                                        if($lignebond[$bi]=~/^\s/){
                                                $extract[0]=$coller[0].$coller[1];
                                                $extract[2]=$extract[1];
                                                $extract[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $extract[0]=$coller[0].$coller[1].$coller[2];
                                                $extract[2]=$extract[1];
                                                $extract[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($lignebond[$bi]=~/^\s/){
                                                $extract[0]=$coller[0];
                                                $extract[2]=$extract[1];
                                                $extract[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $extract[0]=$coller[0].$coller[1].$coller[2];
						$extract[2]=$extract[1];
                                                $extract[1]=$coller[3];
                                        };
                                }
                                elsif(@coller2==4){
                                        $extract[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $extract[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $extract[0]=$coller[0].$coller[1].$coller[2];
                                        $extract[1]=$coller[3].$coller[4].$coller[5];
                                        $extract[2]=$coller[6];
                                };

			$a1=$bi;
			$a2=$extract[0];
			$a3=$extract[1];
			$a4=$extract[2];

			if(($type[$a2] eq 'N.am') && ($type[$a3] eq 'C.2')){
				$a4="am" if($fonc[$a3]=~/ 2-O /);
			};			
			if(($type[$a3] eq 'N.am') && ($type[$a2] eq 'C.2')){
				$a4="am" if($fonc[$a2]=~/ 2-O /);
			};
			
			if(($type[$a2] eq 'O.2') && ($type[$a3] eq 'N.pl3')){
				$a4="2" if($a4=1);
			};
			if(($type[$a2] eq 'N.pl3') && ($type[$a3] eq 'O.2')){
				$a4="2" if($a4=1);
			};
			
			# phosphate
			if(($type[$a2] eq 'P.3') && ($type[$a3] eq 'O.co2')){
				$a4="ar";
			};
			#add25Janv2017
			if(($type[$a2] eq 'O.co2') && ($type[$a3] eq 'P.3')){
				$a4="ar";
			};
	
			# carboxylate
			if(($type[$a2] eq 'C.2') && ($type[$a3] eq 'O.co2')){
				$a4="ar";
			};
			#add25Janv2017
			if(($type[$a2] eq 'O.co2') && ($type[$a3] eq 'C.2')){
                         	$a4="ar";
                        };
			
			# si bond de type 4 dans le sdf alors ar
			if($type[$a3] =~ /.ar/ && $type[$a2] =~ /.ar/ && $fonc[$a3]=~/ 4-/ && $fonc[$a2]=~/ 4-/){
				@det=split(' ',$fonc[$a3]);
				@det1=split(' ',$ifonc[$a3]);
	                        foreach $k (0..@det-1){
	                        	if($det[$k] =~/4-/){
	                        		$a4="ar" if($det1[$k]=$a2);
	                        	};
	                        };
			};
			
			
			if(($type[$a3] =~ /.ar/) && ($type[$a2] =~ /.ar/) ){
				foreach $lc (0..@cyclef-1){
					$tp=' '.$cyclef[$lc].' ';
					$a4="ar" if(($tp=~/ $a2 /) && ($tp=~/ $a3 /) && $cyclefar[$lc] eq "ar");
				};
			};
			if((($type[$a3] =~ /S.2/) && ($type[$a2] =~ /.ar/)) || (($type[$a2] =~ /S.2/) && ($type[$a3] =~ /.ar/))){
				foreach $lc (0..@cyclef-1){
					$tp=' '.$cyclef[$lc].' ';
					$a4="ar" if(($tp=~/ $a2 /) && ($tp=~/ $a3 /)&& $cyclefar[$lc] eq "ar");
				};
			};
			if((($type[$a3] =~ /O.2/) && ($type[$a2] =~ /.ar/)) || (($type[$a2] =~ /O.2/) && ($type[$a3] =~ /.ar/))){
				foreach $lc (0..@cyclef-1){
					$tp=' '.$cyclef[$lc].' ';
					$a4="ar" if(($tp=~/ $a2 /) && ($tp=~/ $a3 /) && $cyclefar[$lc] eq "ar");
				};
			};
			#printf OUTC "%4s%5s%5s $a4\n",$a1,$a2,$a3;
		
			if($framework==0){
				$class_ar++ if($a4=~/ar/);
				$class_1++ if($a4=~/1/);
				$class_2++ if($a4=~/2/);
				$class_3++ if($a4=~/3/);
				$class_am++ if($a4=~/am/);
			}
			elsif($halo!~/ $type[$a3] / && $halo!~/ $type[$a2] / && $exclude_framework!~/ $a3 / && $exclude_framework!~/ $a2 /){
				$class_1++;	
			};	
			
			$typebondmol2[$bi]=$a4;


		};
		#print OUTC "@<TRIPOS>SUBSTRUCTURE\n";
		#print OUTC "      1 mol$moli             1 ****\n";
	#close(OUTC);
	
	
$debug=0;	
	
	## GENERATION KEY

	 %order_sybyl=(
		'C.3',1,
		'C.2',2,
		'C.1',3,
		'C.ar',4,
		'C.cat',5,
		'N.3',6,
		'N.2',7,
		'N.1',8,
		'N.ar',9,
		'N.am',10,
		'N.pl3',11,
		'N.4',12,
		'O.3',13,
		'O.2',14,
		'O.co2',15,
		'O.spc',16,
		'O.t3p',17,
		'S.3',18,
		'S.2',19,
		'S.o',20,
		'S.o2',21,
		'P.3',22,
		'H',23,
		'X',38,
		'H.spc',24,
		'H.t3p',25,
		'F',26,
		'Cl',27,
		'Br',28,
		'I',29,
		'Si',30,
		'LP',31,
		'Du',32,
		'Na',33,
		'K',34,
		'Ca',35,
		'Li',36,
		'Al',37,
		'B',39,
	);

# pour le calcul du centre assymetrique R ou S

	%order_assym=(
		'LP',1,
		'Du',1,
		'H.spc',2,
		'H.t3p',2,
		'H',2,
		'X',25,
		'Li',3,
		'C',4,
		'C.3',4,
		'C.cat',4,
		'C.2',5,
		'C.ar',5,
		'C.1',6,
		'N',7,
		'N.4',7,
		'N.3',7,
		'N.pl3',7,
		'N.am',7,
		'N.2',8,
		'N.ar',8,
		'N.1',9,
		'O',10,
		'O.3',10,
		'O.2',11,
		'O.co2',11,
		'O.spc',10,
		'O.t3p',10,
		'F',12,
		'Na',13,
		'Al',14,
		'Si',15,
		'P',16,
		'P.3',16,
		'S',17,
		'S.3',17,
		'S.2',18,
		'S.o',18,
		'S.o2',19,
		'Cl',20,
		'K',21,
		'Ca',22,
		'Br',23,
		'I',24,
		'B',3,
	);

#print "X = $order_assym{\"X\"}\n";

	@graphe='';
	@node='';
	@no_node='';
	@node_assym='';
	@czore='';
	#print "$istratom atoms\n";
	#print "@type\n";
	#print "@atom\n";
	#print "@atomx\n";

	foreach $bi (1..$istratom){
		        $type[$bi]=$atomx[$bi] if($atomx[$bi] eq 'X' && $eltx eq "X");
			
			if($exclude_framework!~/ $bi /){
				$node[$bi]=" $type[$bi] ";
				#$graphe[$bi]=$type[$bi];
				$no_node[$bi]=" $bi ";
			};
		
		# C assym
			if($type[$bi] eq 'C.3' && $framework==0){
				$node_assym[$bi]=1;
				#print "c no $bi C assymetrique ?\n" if($debug);
			}
			elsif($type[$bi] eq 'P.3' && $framework==0){
				$node_assym[$bi]=0;
				# P.3 with 3 neighbours different from O and 1 =O => asymetric
				if($fonc[$bi] =~/ 2-O / && $fonc[$bi] !~/ 1-O / ){
					$car=$fonc[$bi];
					$nocc=$car=~s/ 2-O / 2-O /g;
					if ($nocc == 1){
						$node_assym[$bi]=1;
						#print "P.3 tetrahydral assymetrique ?\n";
					};
				};
			}
			else{
				$node_assym[$bi]=0;
			};

		# C=C Z or E
			if($type[$bi] eq 'C.2' && $czore[$bi] eq '' && $framework==0){
				@get = split(' ',$ifonc[$bi]);
				@geti = split(' ',$typeb[$bi]);
				$grandzore1="";
				$grandzore2="";
				$nbczore=0;
				$voisinzore1="";
				$voisinzore2="";
				foreach $lc (0..@get-1){
					#works for C=C and C=N-
					if(($type[$get[$lc]] eq 'C.2' && $geti[$lc] == 2) || ($type[$get[$lc]] eq 'N.2' && $geti[$lc] == 2) ){
					#if(($type[$get[$lc]] eq 'C.2' && $geti[$lc] == 2)){
						$nbczore++;
						$voisinzore2=$get[$lc] if($voisinzore1 ne '');
						$voisinzore1=$get[$lc] if($voisinzore1 eq '');
					};
				};
				if($nbczore > 1){ # plusieurs voisins C.2 donc au centre Z et E seront definis par les voisins
					#print "$bi au centre de $voisinzore1 et de $voisinzore2\n";
					$czore[$bi]='';
				}
				elsif($nbczore==1){
					#print "Z ou E ? $bi = $voisinzore1 \n";
					$noboucle=''; # pour empecher de boucler sur des paires deja comparee dans la fonction &compareatom
					################## premier

					$tempifonc=$ifonc[$bi];
					$tempifonc=~s/ $voisinzore1 //g;
					@geti2 = split(' ',$tempifonc);

					$atomtemp1=$atom[$geti2[0]];
					$atomtemp1='X' if($atomx[$geti2[0]] eq 'X' && $eltx eq "X");

					$atomtemp2=$atom[$geti2[1]];
					$atomtemp2='X' if($atomx[$geti2[1]] eq 'X' && $eltx eq "X");

					if($geti2[1] ne ''){

						if($order_assym{"$atomtemp1"} > $order_assym{"$atomtemp2"}){
							$grandzore1=$geti2[0];
						}
						elsif($order_assym{"$atomtemp1"} < $order_assym{"$atomtemp2"}){
							$grandzore1=$geti2[1];
						}
						else{
							#print "$geti2[0] et $geti2[1] => $resultcompar\n";
							$resultcompar=&compareatom($bi,$geti2[0],$bi,$geti2[1]);
							#print "$geti2[0] et $geti2[1] => $resultcompar\n";
							$grandzore1=$geti2[1] if($resultcompar == -1);
							$grandzore1=$geti2[0] if($resultcompar == 1);
						};
						$grandzore1atom=$bi;
					}
					else{
					
					# je l'ai deconnecte car il s'agit en fait d'un axe de chiralite aS ou aR
					$nogo=0;
					if($nogo){
						$searchi=$bi;
						$searchi1=$geti2[0];
						$searchflag=1;
						while($searchflag){

							$tempfonc2=$typeb[$searchi];
							@geti2a = split(' ',$tempfonc2);

							$tempifonc2=$ifonc[$searchi];
							@geti2b = split(' ',$tempifonc2);

							$doyouenter=0;
							foreach $lc (0..@geti2a-1){
								if($geti2b[$lc] ==  $searchi1 && $type[$searchi1] eq 'C.2' && $geti2a[$lc] == 2){
									$doyouenter=1;
									$tempifonc=$ifonc[$searchi1];
									$tempifonc=~s/ $searchi //g;
									@geti2 = split(' ',$tempifonc);
									if($geti2[1] ne ''){
										$searchflag=0;
										$searchi=$searchi1;
									}
									else{
										$searchi=$searchi1;
										$searchi1=$geti2[0];
									};
									last;
								};
							};
							if($doyouenter==0){
								$searchflag=0;
							};
						};

						if($doyouenter==1){
						
							$atomtemp1=$atom[$geti2[0]];
							$atomtemp1='X' if($atomx[$geti2[0]] eq 'X' && $eltx eq "X");

							$atomtemp2=$atom[$geti2[1]];
							$atomtemp2='X' if($atomx[$geti2[1]] eq 'X' && $eltx eq "X");

							if($order_assym{"$atomtemp1"} > $order_assym{"$atomtemp2"}){
								$grandzore1=$geti2[0];
							}
							elsif($order_assym{"$atomtemp1"} < $order_assym{"$atomtemp2"}){
								$grandzore1=$geti2[1];
							}
							else{
								$resultcompar=&compareatom($searchi,$geti2[0],$searchi,$geti2[1]);
								#print "$geti2[0] et $geti2[1] => $resultcompar\n";
								$grandzore1=$geti2[1] if($resultcompar == -1);
								$grandzore1=$geti2[0] if($resultcompar == 1);
							};
							$grandzore1atom=$searchi;
						};
					};
					};

					################## second C.2 or N.2

					$tempifonc=$ifonc[$voisinzore1];
					$tempifonc=~s/ $bi //g;
					@geti2 = split(' ',$tempifonc);

					$atomtemp1=$atom[$geti2[0]];
					$atomtemp1='X' if($atomx[$geti2[0]] eq 'X' && $eltx eq "X");

					$atomtemp2=$atom[$geti2[1]];
					$atomtemp2='X' if($atomx[$geti2[1]] eq 'X' && $eltx eq "X");


					if($geti2[1] ne ''){

						if($order_assym{"$atomtemp1"} > $order_assym{"$atomtemp2"}){
							$grandzore2=$geti2[0];
						}
						elsif($order_assym{"$atomtemp1"} < $order_assym{"$atomtemp2"}){
							$grandzore2=$geti2[1];
						}
						else{
							#print "$geti2[0] et $geti2[1] => $resultcompar\n";
							$resultcompar=&compareatom($voisinzore1,$geti2[0],$voisinzore1,$geti2[1]);
							#print "$geti2[0] et $geti2[1] => $resultcompar\n";
							$grandzore2=$geti2[1] if($resultcompar == -1);
							$grandzore2=$geti2[0] if($resultcompar == 1);
						};
						$grandzore2atom=$voisinzore1;
					}
					elsif($type[$voisinzore1] eq "N.2"){
						#The second neighbour is LP:Lone pair of nitrogen not present in the file
						#Thus, 
						$grandzore2=$geti2[0];
						$grandzore2atom=$voisinzore1;	
					}
					else{
					
					# je l'ai deconnecte car il s'agit en fait d'un axe de chiralite aS ou aR
					$nogo=0;
					if($nogo){
						$searchi=$voisinzore1;
						$searchi1=$geti2[0];
						$searchflag=1;
						while($searchflag){

							$tempfonc2=$typeb[$searchi];
							@geti2a = split(' ',$tempfonc2);

							$tempifonc2=$ifonc[$searchi];
							@geti2b = split(' ',$tempifonc2);

							$doyouenter=0;
							foreach $lc (0..@geti2a-1){
								if($geti2b[$lc] ==  $searchi1 && $type[$searchi1] eq 'C.2' && $geti2a[$lc] == 2){
									$doyouenter=1;
									#print "$searchi1 trouve !\n";
									$tempifonc=$ifonc[$searchi1];
									$tempifonc=~s/ $searchi //g;
									@geti2 = split(' ',$tempifonc);
									if($geti2[1] ne ''){
										$searchflag=0;
										$searchi=$searchi1;
									}
									else{
										$searchi=$searchi1;
										$searchi1=$geti2[0];
									};
									last;
								};
							};
							if($doyouenter==0){
								$searchflag=0;
							};
						};

						if($doyouenter==1){

							$atomtemp1=$atom[$geti2[0]];
							$atomtemp1='X' if($atomx[$geti2[0]] eq 'X' && $eltx eq "X");

							$atomtemp2=$atom[$geti2[1]];
							$atomtemp2='X' if($atomx[$geti2[1]] eq 'X' && $eltx eq "X");

							if($order_assym{"$atomtemp1"} > $order_assym{"$atomtemp2"}){
								$grandzore2=$geti2[0];
							}
							elsif($order_assym{"$atomtemp1"} < $order_assym{"$atomtemp2"}){
								$grandzore2=$geti2[1];
							}
							else{
								#print "$geti2[0] et $geti2[1] => $resultcompar\n";
								$resultcompar=&compareatom($searchi,$geti2[0],$searchi,$geti2[1]);
								#print "$geti2[0] et $geti2[1] => $resultcompar\n";
								$grandzore2=$geti2[1] if($resultcompar == -1);
								$grandzore2=$geti2[0] if($resultcompar == 1);
							};

							$grandzore2atom=$searchi;
							#print "$searchi\n";
						};
					};
					};

					if($grandzore1 ne '' && $grandzore2 ne '' ){
						#print "cote $grandzore1atom grand1=$grandzore1 et cote $grandzore2atom grand2=$grandzore2\n";

						$reponsezoue=&zoue($grandzore1,$grandzore1atom,$grandzore2,$grandzore2atom);

						$czore[$grandzore1atom]=$reponsezoue;
						$czore[$grandzore2atom]=$reponsezoue;

						print "atome $grandzore1atom $czore[$grandzore1atom] et atome $grandzore2atom $czore[$grandzore2atom]\n" if($debug);
					}
					else{
						print "$grandzore1atom et $grandzore2atom ni Z ni E\n" if($debug);
						$czore[$grandzore1atom]="no";
						$czore[$grandzore2atom]="no";
						#print "atome $grandzore1atom $czore[$grandzore1atom] et atome $grandzore2atom $czore[$grandzore2atom]\n";
					};
				};
			};# end C=C Z or E
	};
	
	foreach $bi (1..$istratom){
		if($exclude_framework!~/ $bi /){
		@get = split(' ',$ifonc[$bi]);
		$node_assym[$bi]=$ifonc[$bi] if($node_assym[$bi] == 1);
		foreach $lc (0..@get-1){
			if($exclude_framework!~/ $get[$lc] /){
				if($node[$bi] eq ''){
					$node[$bi]=" $type[$get[$lc]] ";
					$no_node[$bi]=" $get[$lc] ";
				}
				else{
					$node[$bi]=$node[$bi]." $type[$get[$lc]] ";
					$no_node[$bi]=$no_node[$bi]." $get[$lc] ";
				};
			};
		};
		};
	};

	#print "graphe: \n";
	#foreach $bi (1..$istratom){
		#$graphe[$bi]=$type[$bi]."_".$node[$bi];
		#print "$graphe[$bi]\n";
		#print "$bi) $no_node[$bi]\n\n";
	#};
	
	@no_node2=@no_node;   # concatenation feuilles avec blanc
	@node2=@node;
	@no_node3=@no_node; # origin premieres feuilles
	@node3=@node;
	foreach $bi2 (1..$istratom){
		if($exclude_framework!~/ $bi2 /){
		foreach $bi (1..$istratom){
			if($exclude_framework!~/ $bi /){
			@get = split(' ',$no_node2[$bi]);
			$nodetmp='';
			$nodetmp2='';
			foreach $lc (1..@get-1){
				@get2 = split(' ',$node3[$get[$lc]]);
				@get3 = split(' ',$no_node3[$get[$lc]]);
				foreach $lc2 (1..@get2-1){
					if($no_node2[$bi]!~/ $get3[$lc2] / && $nodetmp2!~/ $get3[$lc2] /){
						if($exclude_framework!~/ $get3[$lc2] /){
							$nodetmp=$nodetmp." ".$get2[$lc2]." ";
							$nodetmp2=$nodetmp2." ".$get3[$lc2]." ";
							#print "$lc ($lc2) $get3[$lc2] add\n";
						};
					}
					else{
						#print "$lc ($lc2) $get3[$lc2]\n";
					};
				};
			};
			if($nodetmp ne ''){
				$node2[$bi]=$node2[$bi]." ".$nodetmp;
				$no_node2[$bi]=$no_node2[$bi]." ".$nodetmp2;
				$node[$bi]=$node[$bi]."_".$nodetmp;
				$no_node[$bi]=$no_node[$bi]."_".$nodetmp2;
			};
			};
		};
		};
	};

	#print "graphe: \n";
	foreach $bi (1..$istratom){
		@get=split(' ',$no_node2[$bi]);
		$longget=@get;
		$bip=$bi."_";
		$no_node[$bi]=~s/^ $bi/$bip/;
		#print "carbon no $bi assym ? : $node_assym[$bi]\n" if($node_assym[$bi] != 0);
		#print "$node[$bi]\n";
		#print "$bi) $no_node[$bi] ($longget atoms)\n\n";
	};
	
	foreach $bi (1..$istratom){
		#print "$node[$bi]\n" if($bi==1);
		@get=split('_',$node[$bi]);
		$graphetmp='';
		#rank by order of sybyl atom type

		foreach $bi2 (0..@get-1){
			@get2=split(' ',$get[$bi2]);
			$graphetmp='';
			foreach $bi4 (0..@get2-1){
				$maxorder=100;
				$nomaxorder=-1;
				foreach $bi3 (0..@get2-1){
					if($order_sybyl{"$get2[$bi3]"} < $maxorder && $get2[$bi3] ne ''){
						$maxorder=$order_sybyl{$get2[$bi3]};
						$nomaxorder=$bi3;
					};
				};
				
				$graphetmp=$graphetmp.$get2[$nomaxorder]." " if($nomaxorder != -1);
				$get2[$nomaxorder]='' if($nomaxorder != -1);
				
			};
			#print "$bi $graphetmp\n";
			$graphetmp=~s/ $//;
			$graphetmp=~s/ /-/g;
			$graphe[$bi]=$graphe[$bi]."_".$graphetmp;
			#print "$bi $graphe[$bi]\n";
		};
	};
	
	# compare les 4 feuilles des carbon assymetrique
	$nonassym=0;
	foreach $bi (1..$istratom){
		$flagnosym=0;
		#print "noboucle 1 $noboucle\n";
		@get=split(' ',$node_assym[$bi]);
		if($node_assym[$bi] != 0){
			$longgetifonc=@get;
			if($longgetifonc == 4){
				foreach $i (0 .. @get-2){
					foreach $j ($i+1 .. @get-1){
						if($graphe[$get[$i]] eq $graphe[$get[$j]]){
							$nonassym=1;
							#print "carbon no $bi non assymetric par $get[$i] = $get[$j]\n";
						};
					};
				};
			}
			else{
				$nonassym=1;
				print "C $bi C.3 avec $longgetifonc voisins !!! molecule exclue !\n";
				$weird=1;
				last;
			};
			if($nonassym==1){
				#print "carbon no $bi non assymetric\n";
				$node_assym[$bi]=0;
			}
			else{
				print "CARBON no $bi assymetric ?\n" if($debug);
				$node_assym[$bi]="assym";

				################################################################"
				# rechercher ordre

				$ifonctemp=" ".$ifonc[$bi]." ";
				@getifonc=split(' ',$ifonctemp);

				print "$bi ordonner $ifonctemp\n" if($debug);
				$turn=0;
				$sortiewh=0;
				while(@getifonc >= 1 && $sortiewh==0){
					$petit=$getifonc[0];
					@getifoncold=@getifonc;
					#print "probleme ici @getifonc\n";

					foreach $i (0..@getifonc-1){

						$atomtemp1=$atom[$petit];
						$atomtemp1='X' if($atomx[$petit] eq 'X' && $eltx eq "X");

						$atomtemp2=$atom[$getifonc[$i]];
						$atomtemp2='X' if($atomx[$getifonc[$i]] eq 'X' && $eltx eq "X");

						if($order_assym{"$atomtemp1"} > $order_assym{"$atomtemp2"}){
							$petit=$getifonc[$i];
						};
					};
					
					#print "petit $petit\n";

					$flagidem=0;
					$tocompar="";
					foreach $i (0..@getifonc-1){

						$atomtemp1=$atom[$petit];
						$atomtemp1='X' if($atomx[$petit] eq 'X' && $eltx eq "X");

						$atomtemp2=$atom[$getifonc[$i]];
						$atomtemp2='X' if($atomx[$getifonc[$i]] eq 'X' && $eltx eq "X");

						if($order_assym{"$atomtemp1"} == $order_assym{"$atomtemp2"}){
							$flagidem++;
							$tocompar=$tocompar." $getifonc[$i] ";
						};
					};
					print "petit $petit et tocompar $tocompar\n" if($debug);

					if($flagidem > 1){
						#print "tour $turn : $tocompar a comparer\n";
						$petit="";
						@getordre=split(' ',$tocompar);
						$sortie=0;
						while(@getordre > 1 && $sortie==0){

							foreach $i (0..@getordre-2){
							#print "$i tour $turn : $tocompar a comparer\n";
								foreach $j ($i+1..@getordre-1){
									#print "$j tour $turn : $tocompar a comparer\n";
									$noboucle=''; # pour empecher de boucler sur des paires deja comparee dans la fonction &compareatom
									#print "compareatom $bi $getordre[$i] $bi $getordre[$j]\n";
									$resultcompar=&compareatom($bi,$getordre[$i],$bi,$getordre[$j]);
									#print "noboucle 1 $noboucle\n";
									#print "ok\n";

									print "$getordre[$i] et $getordre[$j] -> $resultcompar\n" if($debug);

									if($resultcompar == 0){
										print "WARNING molecule $nombre_molecule atom $bi assym or not car deux substituants identiques (ambiguities !)\n";
										$sortie=1;
									}
									elsif($resultcompar == -1){ # $getordre[$i] inferieur a $getordre[$j]
										$tocompar=~s/ $getordre[$j] //g;
									}
									else{ # $getordre[$i] superieur a $getordre[$j]
										$tocompar=~s/ $getordre[$i] //g;
									};
								};
								# add 22 sept 2017 (PB boucle infinie sur tocompar empty)
								@getordre=split(' ',$tocompar);
							};
							@getordre=split(' ',$tocompar);
							#print "getordre : @getordre\n";
						};

						$petit=$tocompar;
						$petit=$getordre[0] if($sortie==1);
						$petit=~s/ //g;
						$flagnosym=1 if($sortie==1);
					};

				#	print "tour $turn : le plus petit entre $ifonctemp = $tocompar\n";
					$turn++;
					if($turn==1){
						$atomcd=$petit;
					}
					elsif($turn==2){
						$atomcc=$petit;
					}
					elsif($turn==3){
						$atomcb=$petit;
					}
					elsif($turn==4){
						$atomca=$petit;
					};

					$ifonctemp=~s/ $petit //g;
					@getifonc=split(' ',$ifonctemp);
					
					if(@getifonc eq @getifoncold){
						$sortiewh=1;
						print "Desaccord sur ordre ex: a<b b<c et a>c !\n";
						$node_assym[$bi]=0;
					};
					
					print "$bi ordonner $ifonctemp (petit courant = $petit)\n" if($debug);
				};

				#print "Autour de $bi ordre par type : $atomcd < $atomcc < $atomcb < $atomca\n";
				$order1=" $atomca $atomcb $atomcc ";
				$order2=" $atomcc $atomca $atomcb ";
				$order3=" $atomcb $atomcc $atomca ";

				################################################################
				if($sortiewh==0){
				if($flagnosym==0){

					&assymrous($bi,$ifonc[$bi],$atomcd);
					#print "Autour de ordre clockwise $clockwise\n";

					if($order1 eq $clockwise){
						$node_assym[$bi]="assymR";
					}
					elsif($order2 eq $clockwise){
						$node_assym[$bi]="assymR";
					}
					elsif($order3 eq $clockwise){
						$node_assym[$bi]="assymR";
					}
					else{
						$node_assym[$bi]="assymS";
					};
					print "$clockwise / $order1 ; $order2 ; $order3\n" if($debug);
					print "atom no $bi $node_assym[$bi]\n\n" if($debug);
				}
				else{
					$node_assym[$bi]="tosee";
				};
				};
			};
			$nonassym=0;
		};
	};

	# un second tour sur les atomes assymetriques 'tosee' qui avaient deux substituants identiques
	$nonassym=0;
	foreach $bi (1..$istratom){
		if($node_assym[$bi] eq "tosee" && $weird==0){
			print "re-check atom $bi\n";

			$flagnosym=0;
			#print "noboucle 1 $noboucle\n";

			@get=split(' ',$node_assym[$bi]);


			foreach $i (0 .. @get-2){
				foreach $j ($i+1 .. @get-1){
					if($graphe[$get[$i]] eq $graphe[$get[$j]]){
						$nonassym=1;
						#print "carbon no $bi non assymetric par $get[$i] = $get[$j]\n";
					};
				};
			};

			if($nonassym==1){
				#print "carbon no $bi non assymetric\n";
				$node_assym[$bi]=0;
			}
			else{
				#print "carbon no $bi assymetric ?\n";
				$node_assym[$bi]="assym";

				################################################################"
				# rechercher ordre

				$ifonctemp=" ".$ifonc[$bi]." ";
				@getifonc=split(' ',$ifonctemp);
				#print "ordonner $ifonctemp\n";
				$turn=0;
				while(@getifonc >= 1){
					$petit=$getifonc[0];
					foreach $i (0..@getifonc-1){

						$atomtemp1=$atom[$petit];
						$atomtemp1='X' if($atomx[$petit] eq 'X' && $eltx eq "X");

						$atomtemp2=$atom[$getifonc[$i]];
						$atomtemp2='X' if($atomx[$getifonc[$i]] eq 'X' && $eltx eq "X");

						if($order_assym{"$atomtemp1"} > $order_assym{"$atomtemp2"}){
							$petit=$getifonc[$i];
						};
					};
					$flagidem=0;
					$tocompar="";
					foreach $i (0..@getifonc-1){

						$atomtemp1=$atom[$petit];
						$atomtemp1='X' if($atomx[$petit] eq 'X' && $eltx eq "X");

						$atomtemp2=$atom[$getifonc[$i]];
						$atomtemp2='X' if($atomx[$getifonc[$i]] eq 'X' && $eltx eq "X");

						if($order_assym{"$atomtemp1"} == $order_assym{"$atomtemp2"}){
							$flagidem++;
							$tocompar=$tocompar." $getifonc[$i] ";
						};
					};

					if($flagidem > 1){
						#print "tour $turn : $tocompar a comparer\n";
						$petit="";
						@getordre=split(' ',$tocompar);
						$sortie=0;
						while(@getordre > 1 && $sortie==0){
							foreach $i (0..@getordre-2){
							#print "$i tour $turn : $tocompar a comparer\n";
								foreach $j ($i+1..@getordre-1){
									#print "$j tour $turn : $tocompar a comparer\n";
									$noboucle=''; # pour empecher de boucler sur des paires deja comparee dans la fonction &compareatom
									$resultcompar=&compareatom($bi,$getordre[$i],$bi,$getordre[$j]);
									#print "noboucle 1 $noboucle\n";
									#print "ok\n";
	
									#print "avant tocompar $tocompar \n";
									if($resultcompar == 0){
										print "WARNING molecule $nombre_molecule atom $bi assym or not car deux substituants identiques (ambiguities !)\n";
										#$weird=1;
										$sortie=1;
									}
									elsif($resultcompar == -1){ # $getordre[$i] inferieur a $getordre[$j]
										$tocompar=~s/ $getordre[$j] //g;
									}
									else{ # $getordre[$i] superieur a $getordre[$j]
										$tocompar=~s/ $getordre[$i] //g;
									};
								};
								#print "tocompar $tocompar\n";
								# add 22 sept 2017 (PB boucle infinie sur tocompar empty)
								@getordre=split(' ',$tocompar);
							};
							@getordre=split(' ',$tocompar);
							#print "getordre : @getordre\n";
						};
						$petit=$tocompar;
						$petit=$getordre[0] if($sortie==1);
						$petit=~s/ //g;
						$flagnosym=1 if($sortie==1);
					};

					#print "tour $turn : le plus petit entre $ifonctemp = $tocompar\n";
					$turn++;
					if($turn==1){
						$atomcd=$petit;
					}
					elsif($turn==2){
						$atomcc=$petit;
					}
					elsif($turn==3){
						$atomcb=$petit;
					}
					elsif($turn==4){
						$atomca=$petit;
					};

					$ifonctemp=~s/ $petit //g;
					@getifonc=split(' ',$ifonctemp);
				};

				#print "Autour de $bi ordre par type : $atomcd < $atomcc < $atomcb < $atomca\n";
				$order1=" $atomca $atomcb $atomcc ";
				$order2=" $atomcc $atomca $atomcb ";
				$order3=" $atomcb $atomcc $atomca ";

				################################################################
				if($flagnosym==0){

					&assymrous($bi,$ifonc[$bi],$atomcd);
					#print "Autour de ordre clockwise $clockwise\n";

					if($order1 eq $clockwise){
						$node_assym[$bi]="assymR";
					}
					elsif($order2 eq $clockwise){
						$node_assym[$bi]="assymR";
					}
					elsif($order3 eq $clockwise){
						$node_assym[$bi]="assymR";
					}
					else{
						$node_assym[$bi]="assymS";
					};

					print "atom no $bi $node_assym[$bi]\n\n" if($debug);
				}
				else{
					$node_assym[$bi]=0;
				};
			};
			$nonassym=0;
		};
	};

	# Creation de la cle
	$keyg='';
	#print "graphe: \n";
	$gi=0;
	foreach $bi (1..$istratom){
		@get=split(' ',$no_node2[$bi]);
		$longget=@get;
		$graphe[$bi]=~s/^_//;
		if($framework==0){
			#print "$bi $atom[$bi] $node_assym[$bi]\n$graphe[$bi]\n" if($node_assym[$bi] =~/^assym/ && $debug);
			print "$bi $atom[$bi] $node_assym[$bi]\n" if($node_assym[$bi] =~/^assym/ && $debug);
			$graphe[$bi]=$node_assym[$bi].$graphe[$bi] if($node_assym[$bi] =~/^assym/);
			$graphe[$bi]=$czore[$bi].$graphe[$bi] if($czore[$bi] ne '' && $czore[$bi] ne "no");
			#print "Apres:\n$graphe[$bi]\n" if($node_assym[$bi] =~/^assym/);
		}
		else{
			if($halo=~/ $atom[$bi] /){
				$graphe[$bi]="";
			}
			else{	
			
			$grapheold=$graphe[$bi];
			
			# replace heteroatoms by C and remove H
			$grapheold=~s/\.ar//g;
		        $grapheold=~s/\.2//g;
			$grapheold=~s/\.3//g;
			$grapheold=~s/\.pl3//g;
			$grapheold=~s/\.1//g;
			$grapheold=~s/\.4//g;
			$grapheold=~s/\.cat//g;
			$grapheold=~s/\.am//g;
			$grapheold=~s/\.co2//g;
			$grapheold=~s/\.spc//g;
			$grapheold=~s/\.t3p//g;
			$grapheold=~s/\.o2//g;
			$grapheold=~s/\.o//g;

			$grapheold=~s/N/C/g;
			$grapheold=~s/O/C/g;
			$grapheold=~s/P/C/g;
			$grapheold=~s/S/C/g;
			$grapheold=~s/Si/C/g;
			$grapheold=~s/LP/C/g;
			$grapheold=~s/Du/C/g;
			$grapheold=~s/Na/C/g;
			$grapheold=~s/K/C/g;
			$grapheold=~s/Ca/C/g;
			$grapheold=~s/Li/C/g;
			$grapheold=~s/Al/C/g;
		
			$grapheold=~s/H//g;
			$grapheold=~s/X//g;
			$grapheold=~s/F//g;
			$grapheold=~s/Cl//g;
			$grapheold=~s/Br//g;
			$grapheold=~s/I//g;

			$grapheold=~s/B/C/g;
		
			# remove extra '-', '_' and '*'
			@getframe=split('_',$grapheold);
			$grapheold="";
			foreach $gfi (0..@getframe-1){
				if($getframe[$gfi] ne ""){
					@getframe2=split('-',$getframe[$gfi]);
					$tmpframe="";
					foreach $gfi2 (0..@getframe2-1){
						$tmpframe=$tmpframe."-".$getframe2[$gfi2] if($getframe2[$gfi2] ne "");
					};
					$tmpframe=~s/^-//;
					$tmpframe=~s/-$//;	
					$grapheold=$grapheold."_".$tmpframe if($tmpframe ne "");
				};	
			};
			$grapheold=~s/^_//;
			$grapheold=~s/_$//;
			
			$graphe[$bi]=$grapheold;
			};
		};	
		if($graphe[$bi] ne ""){
			$keyg=$keyg.$graphe[$bi]."*";
			$gi++;
		};
	};
	print "PB in molecule $nombre_molecule: Nb atomes $class_c3 !=  $gi graphes printed in framework mode (excluded atoms: $exclude_framework)\n" if($framework && $class_c3 != $gi);

	#print "$keyg\n";
	$key2=$keyg;
	$key2=~s/\*$//;
	#print "$keyg\n";

	$key1='';

	$key1="$class_c3-$class_c2-$class_c1-$class_car-$class_cat-$class_n3-$class_n2-$class_n1-$class_nar-$class_nam-$class_npl3-$class_n4-$class_o3-$class_o2-$class_oco2-$class_ospc-$class_ot3p-$class_s3-$class_s2-$class_so-$class_so2-$class_p3-$class_h-$class_hspc-$class_ht3p-$class_f-$class_cl-$class_br-$class_i-$class_si-$class_lp-$class_du-$class_na-$class_k-$class_ca-$class_li-$class_al-$class_ar-$class_1-$class_2-$class_3-$class_am";
	
	$ajour=1;

}
elsif($framework){ # acyclic
	$weird=1;
};# if non acyclic


if($weird==0){
	
	#$lkeya=@key_a;
	#print "KEY\n@key_a\n";
	$occurencei=-1;	
     	foreach $i (0..@key_a-1){
     		#print "($lkeya) KEY $i\n";
     		
     		if($key1 eq $key_a[$i]){
     			print "Cle primaire equivalente entre $i $origin[$i] et $fileopen2 $nombre_molecule \n" if($debug);
     		        $sortie=0;
     			@get=split('\*',$key_b[$i]);
     			@get2=split('\*',$key2);
     			$lget=@get;
     			$lget2=@get2;
     			
     			#print "key2\n $key2\n";
     			#print "key1\n $key_b[$i]\n";
     			
     			#print "$key1\n";
     			#print "$key_b[$i]\n";
     			
     			if($lget==$lget2){
     				foreach $j (0..@get-1){
     					$branche=0;
     				
     					foreach $j2 (0..@get2-1){
     						print "$j $j2 Cle secondaire equivalente ? $get[$j] == $get2[$j2]\n" if($debug);
     						if($get[$j] eq $get2[$j2] && $get2[$j2] ne ""){
     						        print "$j $j2 branche equivalente $get[$j] == $get2[$j2]\n" if($debug);
     							$get2[$j2]="";
     							$branche=1;
     							last;
     						};
     					};
     					if($branche==0){
     						print "$j branche $get[$j] non match\n" if($debug);
     						$sortie=1;
     						last;
     					};
     				};
     		        	if($sortie){
     		        		print "Cle secondaire differente cle \n" if($debug);
     		        	}
     		       	 	else{
     		        		print "Cle secondaire equivalente entre $origin[$i] et $fileopen2 $nombre_molecule\n" if($debug);
     		        		$ajour=0;
     		        		$same=1;
     		        		$laquelle=$origin[$i];
					$occurencei=$i;
     		        	};
     		
     		        }
     		        elsif($framework==0){
     		        	print "$nombre_molecule PB : meme cle primaire que $origin[$i] mais pas le meme nombre d'atomes ?\n";
     		        };
     		}
     		else{
     			print "Cle primaire differente cle $origin[$i]\n" if($debug);
     		};
     		
     		if($ajour==0){
     			last;
     		};
     };

#	if(!-e "$filekey" || -z "$filekey"){
#		system("touch $filekey");
#	};	
	
	#open(OUT,">>$filekey");
	#print OUT "$key1*$key2 $fileopen2 $nombre_molecule\n";
	#close(OUT);
     	#print "Mise a jour de $filekey par $fileopen2 $nombre_molecule\n";
	#print "$fileopen2 $nombre_molecule diff\n";
     	
	#$key_b[@key_b]=$key2;
	#$key_a[@key_a]=$key1;
	#$origin[@origin]="$fileopen2 $nombre_molecule";
     	
	#}
	#elsif($ajour){

	$info2=$nombre_molecule;	
	$info2="NO_".$nombre_molecule."_MDLNUMBER_".$info if($info ne "");	
	
	if($ajour && $filekey ne ""){ #$filekey ne "" add for SKETCH.pl that uses &convert	
		open(OUT,">>$filekey");
     		print OUT "$key1*$key2 $fileopen2 $info2\n";
     		close(OUT);
     		#print "Mise a jour de $filekey $fileopen2 $info2\n";
     		print "$fileopen2 $info2 diff\n" if($debug);
     		
		if($key_b[0] eq ""){
			$bi=0;
		}
		else{	
			$bi=@key_b;
		};	
		
     		$key_b[$bi]=$key2;
     		$key_a[$bi]=$key1;
		$occurence[$bi]=1;
     		$origin[$bi]="$fileopen2 $info2";
     	}
     	else{
		$samei=$occurencei;
     		print "$fileopen2 $info2 same $laquelle key tab_indice_no $occurencei\n" if($filekey ne "");# if($filekey ne "") add for SKETCH.pl that uses &convert	
		$occurence[$occurencei]=$occurence[$occurencei]+1 if($occurencei > -1);
     	};
	
};#if $weird==0;

	# used by MATCH_KEY.pl
	$wholekey="$key1*$key2";
	$wholekey;
};



######################################################################
######################################################################



sub cyclesdf2{

	&cyclesdf;

	# check aromaticity
	
       	@cyclefar='';
	foreach $lc (0..@cyclef-1){       
		@arom=split(' ',$cyclef[$lc]);
		#print "test aromaticity of $cyclef[$lc]\n";
		if(@arom == 5 || @arom == 6){
			$longcycle=@arom;
			$nbatomarpartial=0;
			foreach $lc2 (0..@arom-1){
				$narp=" ".$fonc[$arom[$lc2]]." ";
				$nocc1=$narp=~s/1-/1-/g;
				$nocc2=$narp=~s/2-/2-/g;
				if($atom[$arom[$lc2]] eq 'C'){
					$nbatomarpartial++ if($nocc1==2 && $nocc2==1);
				}
				elsif($atom[$arom[$lc2]] eq 'O'){	
					$nbatomarpartial++ if($nocc1==2 && @arom==5);
				}
				elsif($atom[$arom[$lc2]] eq 'N'){
					if($nocc1==1 && $nocc2==1){
						$nbatomarpartial++;
					}
					elsif($nocc1==3 && @arom==5){
						$nbatomarpartial++;
					};
				}
				elsif($atom[$arom[$lc2]] eq 'S'){
					$nbatomarpartial++ if($nocc1==2 && @arom==5);	
				};
			};
			if($nbatomarpartial==$longcycle){
				$cyclefar[$lc]="ar";
				#print "$cyclef[$lc] aromatique\n";
				foreach $lc2 (0..@arom-1){
					$type[$arom[$lc2]]="N.ar" if($atom[$arom[$lc2]] eq 'N' && $type[$arom[$lc2]] ne 'N.am');
					$type[$arom[$lc2]]="C.ar" if($atom[$arom[$lc2]] eq 'C');
					$type[$arom[$lc2]]="S.2" if($atom[$arom[$lc2]] eq 'S');
					$type[$arom[$lc2]]="O.2" if($atom[$arom[$lc2]] eq 'O');
				};
			};	
		};
	};		
};


#################################################################################################################################


sub assymrous{
	local($cassym,$abcd,$petit)=@_;

	@getcassym=split(' ',$abcd);
	$cassyma=$getcassym[0];
	$cassymb=$getcassym[1];
	$cassymc=$getcassym[2];
	$cassymd=$getcassym[3];

	#print "sub assymrous : $cassyma $cassymb $cassymc $cassymd et petit = $petit\n";

	# translation origine du carbone : 0.0
	$cassymx=$strx[$cassym]-$strx[$cassym];
	$cassymy=$stry[$cassym]-$stry[$cassym];
	$cassymz=$strz[$cassym]-$strz[$cassym];

	# translation de a
	$cassymax=$strx[$cassyma]-$strx[$cassym];
	$cassymay=$stry[$cassyma]-$stry[$cassym];
	$cassymaz=$strz[$cassyma]-$strz[$cassym];

	# translation de b
	$cassymbx=$strx[$cassymb]-$strx[$cassym];
	$cassymby=$stry[$cassymb]-$stry[$cassym];
	$cassymbz=$strz[$cassymb]-$strz[$cassym];

	# translation de c
	$cassymcx=$strx[$cassymc]-$strx[$cassym];
	$cassymcy=$stry[$cassymc]-$stry[$cassym];
	$cassymcz=$strz[$cassymc]-$strz[$cassym];

	# translation de d
	$cassymdx=$strx[$cassymd]-$strx[$cassym];
	$cassymdy=$stry[$cassymd]-$stry[$cassym];
	$cassymdz=$strz[$cassymd]-$strz[$cassym];

	if($petit == $getcassym[0]){
		$petitx=$cassymax;
		$petity=$cassymay;
		$petitz=$cassymaz;
	}
	elsif($petit == $getcassym[1]){
		$petitx=$cassymbx;
		$petity=$cassymby;
		$petitz=$cassymbz;
	}
	elsif($petit == $getcassym[2]){
		$petitx=$cassymcx;
		$petity=$cassymcy;
		$petitz=$cassymcz;
	}
	elsif($petit == $getcassym[3]){
		$petitx=$cassymdx;
		$petity=$cassymdy;
		$petitz=$cassymdz;
	};

 #Vecteurs
	$ux=$cassymx-$petitx;
	$uy=$cassymy-$petity;
	$uz=$cassymz-$petitz;
	$nu=sqrt($ux**2+$uy**2+$uz**2);

	$vx=$cassymx-0;
	$vy=$cassymy-0;
	$vz=$cassymz+2;
	$nv=sqrt($vx**2+$vy**2+$vz**2);

# Angle entre les deux vecteurs a fusionner
	$cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);
	$testcos=1-$cosangle*$cosangle;
	$testcos=sqrt($testcos*$testcos);
	$teta=atan2(sqrt($testcos),$cosangle);
	$tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;
	#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz)\n";
	#print "AVANT angle(vecteur1,vecteur2) $tetab degres\n";

# Vecteur normes
	$ux=$ux/$nu;
	$uy=$uy/$nu;
	$uz=$uz/$nu;

	$vx=$vx/$nv;
	$vy=$vy/$nv;
	$vz=$vz/$nv;

# Calcul des vecteurs orthogonaux
	# W :
	$wx=($uy*$vz)-($uz*$vy);
	$wy=($uz*$vx)-($ux*$vz);
	$wz=($ux*$vy)-($uy*$vx);
	$nw=sqrt($wx**2+$wy**2+$wz**2);


if($nw != 0){

	$wx=$wx/$nw;
	$wy=$wy/$nw;
	$wz=$wz/$nw;

	# K :
	$kx=($wy*$uz)-($wz*$uy);
	$ky=($wz*$ux)-($wx*$uz);
	$kz=($wx*$uy)-($wy*$ux);
	$nk=sqrt($kx**2+$ky**2+$kz**2);

	$kx=$kx/$nk;
	$ky=$ky/$nk;
	$kz=$kz/$nk;

	########################################
	# ROTATION
	# valeur de l'angle

	$cosangle=(($vx*$kx)+($vy*$ky)+($vz*$kz))/($nv * $nk);
#	$angle=atan2(sqrt(1-$cosangle*$cosangle),$cosangle)*180/3.14159;
#	print"angle(v,k) $angle (cos=$cosangle)\n";
	$teta=sqrt($teta*$teta);
	$teta=-1*$teta if($cosangle < 0);


	#######################################
	# pour toutes les coordonnees du vecteur 1:
	# On a besoin de teta, u, k et w

		# a
		$ncorx=$cassymax;
		$ncory=$cassymay;
		$ncorz=$cassymaz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymax=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymay=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymaz=($xrot*$uz + $yrot*$kz + $zrot*$wz);

		# b
		$ncorx=$cassymbx;
		$ncory=$cassymby;
		$ncorz=$cassymbz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymbx=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymby=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymbz=($xrot*$uz + $yrot*$kz + $zrot*$wz);

		# c
		$ncorx=$cassymcx;
		$ncory=$cassymcy;
		$ncorz=$cassymcz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymcx=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymcy=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymcz=($xrot*$uz + $yrot*$kz + $zrot*$wz);

		# d
		$ncorx=$cassymdx;
		$ncory=$cassymdy;
		$ncorz=$cassymdz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymdx=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymdy=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymdz=($xrot*$uz + $yrot*$kz + $zrot*$wz);
}
else{

		# inverse les coordonnees
		$cassymax=$cassymax*-1;
		$cassymay=$cassymay*-1;
		$cassymaz=$cassymaz*-1;

		$cassymbx=$cassymbx*-1;
		$cassymby=$cassymby*-1;
		$cassymbz=$cassymbz*-1;

		$cassymcx=$cassymcx*-1;
		$cassymcy=$cassymcy*-1;
		$cassymcz=$cassymcz*-1;

		$cassymdx=$cassymdx*-1;
		$cassymdy=$cassymdy*-1;
		$cassymdz=$cassymdz*-1;


};

## VERIFIE ANGLE

	if($petit == $getcassym[0]){
		$petitx=$cassymax;
		$petity=$cassymay;
		$petitz=$cassymaz;
	}
	elsif($petit == $getcassym[1]){
		$petitx=$cassymbx;
		$petity=$cassymby;
		$petitz=$cassymbz;
	}
	elsif($petit == $getcassym[2]){
		$petitx=$cassymcx;
		$petity=$cassymcy;
		$petitz=$cassymcz;
	}
	elsif($petit == $getcassym[3]){
		$petitx=$cassymdx;
		$petity=$cassymdy;
		$petitz=$cassymdz;
	};

	$ux=$cassymx-$petitx;
	$uy=$cassymy-$petity;
	$uz=$cassymz-$petitz;
	$nu=sqrt($ux**2+$uy**2+$uz**2);

	$vx=$cassymx-0;
	$vy=$cassymy-0;
	$vz=$cassymz-2;
	$nv=sqrt($vx**2+$vy**2+$vz**2);

	# Angle entre les deux vecteurs a fusionner
	$cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);
	$testcos=1-$cosangle*$cosangle;
	$testcos=sqrt($testcos*$testcos);
	$teta=atan2(sqrt($testcos),$cosangle);
	$tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;
	#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz)\n";
	#print "APRES angle(vecteur1,vecteur2) $tetab degres\n";
## FIN VERIFIE ANGLE

 	#print "$cassymax $cassymay $cassymaz\n";
 	#print "$cassymbx $cassymby $cassymbz\n";
 	#print "$cassymcx $cassymcy $cassymcz\n";
 	#print "$cassymdx $cassymdy $cassymdz\n";
	
	@tabcassymx='';
	@tabcassymy='';
	@tabcassymz='';
	
	$tabcassymx[0]=$cassymax;
	$tabcassymy[0]=$cassymay;
	$tabcassymz[0]=$cassymaz;

	$tabcassymx[1]=$cassymbx;
	$tabcassymy[1]=$cassymby;
	$tabcassymz[1]=$cassymbz;

	$tabcassymx[2]=$cassymcx;
	$tabcassymy[2]=$cassymcy;
	$tabcassymz[2]=$cassymcz;

	$tabcassymx[3]=$cassymdx;
	$tabcassymy[3]=$cassymdy;
	$tabcassymz[3]=$cassymdz;


	# cherche le nombre de x > 0
	$nbxplus=0;
	$nbxplus1='';
	$nbxplus2='';
	foreach $rsi (0..@tabcassymx-1){
		if($tabcassymx[$rsi] > 0 && $petit != $getcassym[$rsi]){
			$nbxplus++;
			$nbxplus2=$rsi if($nbxplus2 eq '' && $nbxplus1 ne '');
			$nbxplus1=$rsi if($nbxplus1 eq '');
		};
	};
	print  "WARNING molecule no $nombre_molecule C_ASSYM no $cassym has three x coordinates > 0 ? molecule exclue !\n" if($nbxplus > 2);
	$weird=1 if($nbxplus > 2);

	if($nbxplus > 1){
		if($tabcassymy[$nbxplus1] > $tabcassymy[$nbxplus2]){
			$plusplus=$getcassym[$nbxplus1];
			$plusmoins=$getcassym[$nbxplus2];
		}
		else{
			$plusplus=$getcassym[$nbxplus2];
			$plusmoins=$getcassym[$nbxplus1];
		};
	}
	else{
		$plusplus=$getcassym[$nbxplus1];

		# cherche y le plus petit
		$smally=1000;
		$smallyi='';
		foreach $rsi (0..@tabcassymy-1){
		if($tabcassymy[$rsi] < $smally && $petit != $getcassym[$rsi] && $plusplus != $getcassym[$rsi]){
			$smally=$tabcassymy[$rsi];
			$smallyi=$rsi;
		};
		};
		$plusmoins=$getcassym[$smallyi];
	};


	# cherche le troisieme !
	$moinsplus=$getcassym[0] if($plusplus != $getcassym[0] && $plusmoins != $getcassym[0] && $petit != $getcassym[0]);
	$moinsplus=$getcassym[1] if($plusplus != $getcassym[1] && $plusmoins != $getcassym[1] && $petit != $getcassym[1]);
	$moinsplus=$getcassym[2] if($plusplus != $getcassym[2] && $plusmoins != $getcassym[2] && $petit != $getcassym[2]);
	$moinsplus=$getcassym[3] if($plusplus != $getcassym[3] && $plusmoins != $getcassym[3] && $petit != $getcassym[3]);

	#print "++=$plusplus +-=$plusmoins -+=$moinsplus\n";

	$clockwise=" $plusplus $plusmoins $moinsplus ";
};




###############################################################################################
###############################################################################################


sub compareatom{
	# Feb 2017 be carefull ! only variables in my() are conserved during recursivity !!
	my($ati,$ai,$atj,$aj,$listegrandi,$listegrandj,$dejacompareri,$dejacomparerj,$mmjkept,$mmikept,$rescompar)=@_;#modif

	$debog=0;
	print "ENTREZ compare $ai (ati = $ati) et $aj (atj = $atj) [noboucle en entrant = $noboucle]\n" if($debog);
	print "\t ai=$ai ($ifonc[$ai]) ; aj=$aj ($ifonc[$aj])\n" if($debog);

if($ai != $aj){
	$dejacompareri=" ".$ifonc[$ai]." ";#modif
	$dejacompareri=~s/ $ati //;
	$dejacomparerib=$dejacompareri;#modif
	$dejacomparerib=~s/ //g;#modif
	
	$dejacomparerj=" ".$ifonc[$aj]." ";#modif
	$dejacomparerj=~s/ $atj //;
	$dejacomparerjb=$dejacomparerj;#modif
	$dejacomparerjb=~s/ //g;#modif
	$rescompar=0;#modif

	$plusati=0;
	$plusatj=0;

	while($rescompar == 0 && $dejacomparerib ne '' && $dejacomparerjb ne ''){#modif
	
		print "\t while $ai : $dejacompareri / $aj : $dejacomparerj\n" if($debog);

		#print "nobouble $noboucle\n";

		# Cherche l'element le plus fort de chaque cote puis compare
		@cpfi=split(' ',$typeb[$ai]);
		@cpi=split(' ',$ifonc[$ai]);
		$sumai=0;
		#print "\t \t $ai type($type[$ai]) typeb($typeb[$ai]) ifonc($ifonc[$ai]) ati = $ati plusati = $plusati\n";
		foreach $ic (0..@cpi-1){

			$atomtemp1=$atom[$cpi[$ic]];
			$atomtemp1='X' if($atomx[$cpi[$ic]] eq 'X' && $eltx eq "X");
			#print "\t\t\tlook on atom $cpi[$ic] atom=$atomtemp1 sumai=$sumai dejacompareri=$dejacompareri\n";

			if($cpi[$ic] != $ati && $dejacompareri =~ / $cpi[$ic] /){ #modif
				#print "\t\t\t enter sumai=$sumai\n";
				if($order_assym{"$atomtemp1"} > $sumai){
					$sumai=$order_assym{"$atomtemp1"};
				};
				#print "\t\t\t after sumai=$sumai\n";
			}
			elsif($cpi[$ic] == $ati && $plusati==0){
				if($order_assym{"$atomtemp1"} > $sumai){
					if($cpfi[$ic] == 2 || $cpfi[$ic] == 3){
						$sumai=$order_assym{"$atomtemp1"};
						$plusati=$sumai;
					};
				};
			};

		};

		#print "\t \t sumai $sumai and plusati=$plusati for $ati\n";

		$listegrandi="";
		$nbmuli=0;
		foreach $ic (0..@cpi-1){
		
			$atomtemp1=$atom[$cpi[$ic]];
			$atomtemp1='X' if($atomx[$cpi[$ic]] eq 'X' && $eltx eq "X");

			if($cpi[$ic] != $ati && $dejacompareri =~ / $cpi[$ic] /){#modif
				if($order_assym{"$atomtemp1"} == $sumai){
					$listegrandi=$listegrandi." $cpi[$ic] ";
					$nbmuli=$nbmuli+1 if($cpfi[$ic] == 2);
					$nbmuli=$nbmuli+2 if($cpfi[$ic] == 3);
				};
			}
			elsif($cpi[$ic] == $ati){
				if($order_assym{"$atomtemp1"} == $sumai){
					if($cpfi[$ic] == 2 || $cpfi[$ic] == 3){
						$nbmuli=$nbmuli+1 if($cpfi[$ic] == 2);
						$nbmuli=$nbmuli+2 if($cpfi[$ic] == 3);
					};
				};
			};

		};

		print "\t \t listgrandi $listegrandi\n" if($debog);

		@cpfj=split(' ',$typeb[$aj]);
		@cpj=split(' ',$ifonc[$aj]);
		$sumaj=0;
		#print "\t \t $aj type($type[$ai]) typeb($typeb[$aj]) ifonc($ifonc[$aj]) atj = $atj plusatj = $plusatj\n";

		foreach $ic (0..@cpj-1){

			$atomtemp1=$atom[$cpj[$ic]];
			$atomtemp1='X' if($atomx[$cpj[$ic]] eq 'X' && $eltx eq "X");
			#print "\t\t\tlook on atom $cpj[$ic] type=$atomtemp1 sumaj=$sumaj dejacomparerj=$dejacomparerj\n";

			if($cpj[$ic] != $atj && $dejacomparerj =~ / $cpj[$ic] /){#modif
				#print "\t\t\t enter sumaj=$sumaj\n";
				if($order_assym{"$atomtemp1"} > $sumaj){
					$sumaj=$order_assym{"$atomtemp1"};
				};
				#print "\t\t\t after sumaj=$sumaj\n";
			}
			elsif($cpj[$ic] == $atj && $plusatj==0){
				#print "\t\t\t $cpj[$ic] = $atj sumaj=$sumaj\n";
				if($order_assym{"$atomtemp1"} > $sumaj){
					#print "\t\t\tenter cpfj = $cpfj[$ic] \n";
					if($cpfj[$ic] == 2 || $cpfj[$ic] == 3){
						$sumaj=$order_assym{"$atomtemp1"};
						$plusatj=$sumaj;
                                        };

				};
 				#print "\t\t\t out sumaj=$sumaj plusatj=$plusatj\n";
			};

		};
		#print "\t \t sumaj $sumaj and plusatj=$plusatj for $atj\n";

		$listegrandj="";
		$nbmulj=0;
		foreach $ic (0..@cpj-1){

			$atomtemp1=$atom[$cpj[$ic]];
			$atomtemp1='X' if($atomx[$cpj[$ic]] eq 'X' && $eltx eq "X");


			if($cpj[$ic] != $atj && $dejacomparerj =~ / $cpj[$ic] /){
				if($order_assym{"$atomtemp1"} == $sumaj){
					$listegrandj=$listegrandj." $cpj[$ic] ";
					$nbmulj=$nbmulj+1 if($cpfj[$ic] == 2);
					$nbmulj=$nbmulj+2 if($cpfj[$ic] == 3);
				};
			}
			elsif($cpj[$ic] == $atj){
				if($order_assym{"$atomtemp1"} == $sumaj){
					if($cpfj[$ic] == 2 || $cpfj[$ic] == 3){
						$nbmulj=$nbmulj+1 if($cpfj[$ic] == 2);
						$nbmulj=$nbmulj+2 if($cpfj[$ic] == 3);
					};
				};
			};

		};
		print "\t \t listgrandj $listegrandj\n" if($debog);

		@mmir=split(' ',$listegrandi);
		$nbmmir=@mmir;
		@mmjr=split(' ',$listegrandj);
		$nbmmjr=@mmjr;

		print "\tcompare $ai et $aj => leur plus gd substituant de valeur respective sumai = $sumai ($nbmmir substituants + $nbmuli par ajout des liaisons multiples) et sumaj = $sumaj ($nbmmjr substituants + $nbmulj par ajout des liaisons multiples)\n" if($debog);

		#print "listb $typeb[$ai] et $typeb[$aj] \n";
		#print "compare $ai et $aj => leurs plus gd substituant $listegrandi /  $listegrandj\n";
		
		$rescompar=0;
		if($sumaj > $sumai){
			$rescompar=-1;
			print "$ai < $aj (plus fort)\n" if($debog);
		}
		elsif($sumai > $sumaj){
			$rescompar=1;
			print "$ai (plus fort) > $aj \n" if($debog);
		}
		else{

			# meme type d'element donc cherche le plus grand nombre de cet element
			@mmi=split(' ',$listegrandi);
			$nbmmi=@mmi+$nbmuli;

			@mmj=split(' ',$listegrandj);
			$nbmmj=@mmj+$nbmulj;

			if($nbmmi < $nbmmj){
				$rescompar=-1;
				print "$ai < $aj (en plus grand nombre)\n" if($debog);
			}
			elsif($nbmmi > $nbmmj){
				$rescompar=1;
				print "$ai (en plus grand nombre) > $aj \n" if($debog);
			}
			else{
				print "Meme nombre : Compare feuille listgrandi $listegrandi et listgrandj $listegrandj\n" if($debog);
				if($mmi[0] ne '' && $mmj[0] ne ''){
					$mmjkept=join(' ',@mmj);
					$mmikept=join(' ',@mmi);

					# la plus grande feuille i
					#print "\t le plus gd substituant porte par $ai parmi ($listegrandi)  \n";

					while(@mmi > 1){
						#print "\t liste i ($listegrandi)\n";

						$atom1et2=$mmi[0]."-".$mmi[1];
						if($noboucle !~ / $atom1et2 /){
							#print "noboucle1 $noboucle et $atom1et2\n";
							$noboucle=$noboucle." ".$mmi[0]."-".$mmi[1]." ".$mmi[1]."-".$mmi[0]." ";

							$mmikept=join(' ',@mmi);
							print "resfeuillei $ai $mmi[0] $ai $mmi[1]\n" if($debog);
							$resfeuille=&compareatom($ai,$mmi[0],$ai,$mmi[1]);
							#print "i compar $mmi[0] $mmi[1] $rescompar\n";
							@mmi=split(' ',$mmikept);

							if($resfeuille==-1){
								$listegrandi=~s/ $mmi[0] //g;
							}
							elsif($resfeuille==1){
								$listegrandi=~s/ $mmi[1] //g;
							}
							else{
								$listegrandi=~s/ $mmi[1] //g;
							};
							@mmi=split(' ',$listegrandi);
						}
						else{
							print "i mmi mis a nul \n" if($debog);
							@mmi='';
						};
					};
					$grandi=$mmi[0];
					print "\t I) rescompar=$rescompar le plus gd substituant porte par $ai parmi ($listegrandi) est $grandi ($listegrandi)\n" if($debog);

					# la plus grande feuille j
					#print "\t le plus gd substituant porte par $aj parmi ($listegrandj)  \n";

					@mmj=split(' ',$mmjkept);
					while(@mmj > 1){
						#print "\t liste j ($listegrandj)\n" if($debog);
						print "\t liste j ($listegrandj) et boucle $noboucle\n" if($debog);
						$atom1et2=$mmj[0]."-".$mmj[1];
						if($noboucle !~/ $atom1et2 /){

							$noboucle=$noboucle." ".$mmj[0]."-".$mmj[1]." ".$mmj[1]."-".$mmj[0]." ";
							#print "noboucle2 $noboucle et $atom1et2\n";
							#print "A $listegrandj et @mmj\n";

							$mmjkept=join(' ',@mmj);
							print "resfeuillej $aj $mmj[0] $aj $mmj[1]\n" if($debog);
							$resfeuille=&compareatom($aj,$mmj[0],$aj,$mmj[1]);
							@mmj=split(' ',$mmjkept);

							#print "j compar $mmj[0] $mmj[1] $rescompar\n";
							#print "res $resfeuille\n";

							if($resfeuille==-1){
								$listegrandj=~s/ $mmj[0] //g;
							}
							elsif($resfeuille==1){
								$listegrandj=~s/ $mmj[1] //g;
							}
							else{
								$listegrandj=~s/ $mmj[1] //g;
							};

							@mmj=split(' ',$listegrandj);
							#print "B $listegrandj et @mmj\n";
						}
						else{
							print "j mmj mis a nul \n" if($debog);
							@mmj='';
						};
					};
					$grandj=$mmj[0];
					print "\t J) rescompar=$rescompar le plus gd substituant porte par $aj parmi ($listegrandj) est $grandj ($listegrandj)\n" if($debog);

					# Compare i et j
					
					print "check recursivity is it good ? : $ai $listegrandi ; $aj $listegrandj\n" if($debog);
					#errors were observed so recalculate beacause not in my() thus not conserved
						@mmi=split(' ',$listegrandi);
						$grandi=$mmi[0];
						@mmj=split(' ',$listegrandj);
						$grandj=$mmj[0];

					if($grandi ne '' && $grandj ne ''){
						$atom1et2=$grandi."-".$grandj;
						if($noboucle !~/ $atom1et2 /){
							#print "noboucle3 $noboucle et $atom1et2\n";
							$noboucle=$noboucle." ".$grandi."-".$grandj." ".$grandj."-".$grandi." ";
							print "resfeuille3 $ai $grandi $aj $grandj\n" if($debog);
							$rescompar=&compareatom($ai,$grandi,$aj,$grandj);
							#print "compar final compar $grandi $grandj $rescompar\n";
						};
					#}
					#else{
					#	print "grandi et/ou grandj null car deja vu(s) !\n";
					};
					print "resultat entre $ai et $aj : $rescompar\n" if($debug);
					#print "RETOUR a $ai et $aj\n";
				}
				elsif($mmi[0] ne ''){
					$rescompar=1;
					#print "resultat $ai le plus grand car $aj pas de liste associee\n";
				}
				elsif($mmj[0] ne ''){
					$rescompar=-1;
					#print "resultat $aj le plus grand car $ai pas de liste associee\n";
				}
				else{
					#print "mmi et mmj null\n";
					#$rescompar=0;
					if($plusatj != $sumaj){
						$plusatj=0;
					};

					if($plusati != $sumai){
						$plusati=0;
					};
				};
			};
		};#if sumai=sumaj

		#print "\t liste i ($listegrandi)\n";
		#print "\t liste j ($listegrandj)\n";
		#print "SORTIE\n";

		if($rescompar == 0){#modif
			#print "bien entrer $dejacompareri / $dejacomparerj et list $listegrandi / $listegrandj\n";
			@mmi=split(' ',$listegrandi);#modif
			foreach $ic (0..@mmi-1){#modif
				$dejacompareri=~s/ $mmi[$ic] //;#modif
			};#modif
			$dejacomparerib=$dejacompareri;#modif
			$dejacomparerib=~s/ //g;#modif

			@mmj=split(' ',$listegrandj);#modif
			foreach $ic (0..@mmj-1){#modif
				$dejacomparerj=~s/ $mmj[$ic] //;#modif
			};#modif
			$dejacomparerjb=$dejacomparerj;#modif
			$dejacomparerjb=~s/ //g;#modif
			#print"bien sortie ($ai/$aj) i : $dejacompareri / j : $dejacomparerj\n";
		};
		print "rescompar $rescompar\n" if($debug);

	};#modif end of while


	if($rescompar == 0 && (($dejacomparerib eq '' && $dejacomparerjb ne '') || ($dejacomparerib ne '' && $dejacomparerjb eq ''))){
		print "probleme de nombre de substituants de $ai [$dejacomparerib] ou $aj [$dejacomparerjb] (rescompar=$rescompar) la molecule a quelque chose de bizarre !!\n" if($debog);
		#$weird=1;
		if($dejacomparerib eq '' && $dejacomparerjb ne ''){
			$rescompar =-1;
		}
		elsif($dejacomparerib ne '' && $dejacomparerjb eq ''){
			$rescompar =1;
		};
	};

	if($rescompar == -1){
		print "resultat $ai < $aj \n" if($debog);
	}
	elsif($rescompar == 1){
		print "resultat $ai > $aj \n" if($debog);
	}
	else{
		print "$node_assym[$ai] et $node_assym[$aj]\n" if($debog);
		if($node_assym[$ai] eq "assymS" && $node_assym[$aj] eq "assymR"){
			$rescompar =-1;
		}
		elsif($node_assym[$aj] eq "assymS" && $node_assym[$ai] eq "assymR"){
			$rescompar =1;
		};

		if($rescompar == -1){
			print "resultat $ai < $aj car s < r\n" if($debog);
		}
		elsif($rescompar == 1){
			print "resultat $ai > $aj  car r > s \n" if($debog);
		}
		else{
			print "resultat $ai = $aj \n" if($debog);
		};
	};
}
else{
	#print "compare les meme atomes !\n";
	$rescompar =0;
};

	print "SORTIE $ai / $aj\n" if($debog);
	return $rescompar;


};

###############################################################################################
###############################################################################################


sub zoue{

	local($cassyma,$cassymb,$cassymc,$cassymd)=@_;

	#print "sub Z ou E : $cassyma $cassymb $cassymc $cassymd \n";

	$reseouz=0;

	# translation de a
	$cassymax=$strx[$cassyma]-$strx[$cassymb];
	$cassymay=$stry[$cassyma]-$stry[$cassymb];
	$cassymaz=$strz[$cassyma]-$strz[$cassymb];


	# translation de c
	$cassymcx=$strx[$cassymc]-$strx[$cassymb];
	$cassymcy=$stry[$cassymc]-$stry[$cassymb];
	$cassymcz=$strz[$cassymc]-$strz[$cassymb];

	# translation de d
	$cassymdx=$strx[$cassymd]-$strx[$cassymb];
	$cassymdy=$stry[$cassymd]-$stry[$cassymb];
	$cassymdz=$strz[$cassymd]-$strz[$cassymb];


	# translation de b
	$cassymbx=$strx[$cassymb]-$strx[$cassymb];
	$cassymby=$stry[$cassymb]-$stry[$cassymb];
	$cassymbz=$strz[$cassymb]-$strz[$cassymb];


 #Vecteurs
	$ux=$cassymbx-$cassymax;
	$uy=$cassymby-$cassymay;
	$uz=$cassymbz-$cassymaz;
	$nu=sqrt($ux**2+$uy**2+$uz**2);

	# vecteur v en -1 1 0
	$vx=$cassymbx+1;
	$vy=$cassymby-1;
	$vz=$cassymbz-0;
	$nv=sqrt($vx**2+$vy**2+$vz**2);

# Angle entre les deux vecteurs a fusionner
	$cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);
	$testcos=1-$cosangle*$cosangle;
	$testcos=sqrt($testcos*$testcos);
	$teta=atan2(sqrt($testcos),$cosangle);
	$tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;
	#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz)\n";
	#print "AVANT angle(vecteur1,vecteur2) $tetab degres\n";

# Vecteur normes
	$ux=$ux/$nu;
	$uy=$uy/$nu;
	$uz=$uz/$nu;

	$vx=$vx/$nv;
	$vy=$vy/$nv;
	$vz=$vz/$nv;

# Calcul des vecteurs orthogonaux
	# W :
	$wx=($uy*$vz)-($uz*$vy);
	$wy=($uz*$vx)-($ux*$vz);
	$wz=($ux*$vy)-($uy*$vx);
	$nw=sqrt($wx**2+$wy**2+$wz**2);


if($nw != 0){

	$wx=$wx/$nw;
	$wy=$wy/$nw;
	$wz=$wz/$nw;

	# K :
	$kx=($wy*$uz)-($wz*$uy);
	$ky=($wz*$ux)-($wx*$uz);
	$kz=($wx*$uy)-($wy*$ux);
	$nk=sqrt($kx**2+$ky**2+$kz**2);

	$kx=$kx/$nk;
	$ky=$ky/$nk;
	$kz=$kz/$nk;

	########################################
	# ROTATION
	# valeur de l'angle

	$cosangle=(($vx*$kx)+($vy*$ky)+($vz*$kz))/($nv * $nk);
#	$angle=atan2(sqrt(1-$cosangle*$cosangle),$cosangle)*180/3.14159;
#	print"angle(v,k) $angle (cos=$cosangle)\n";
	$teta=sqrt($teta*$teta);
	$teta=-1*$teta if($cosangle < 0);


	#######################################
	# pour toutes les coordonnees du vecteur 1:
	# On a besoin de teta, u, k et w

		# a
		$ncorx=$cassymax;
		$ncory=$cassymay;
		$ncorz=$cassymaz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymax=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymay=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymaz=($xrot*$uz + $yrot*$kz + $zrot*$wz);


		# c
		$ncorx=$cassymcx;
		$ncory=$cassymcy;
		$ncorz=$cassymcz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymcx=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymcy=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymcz=($xrot*$uz + $yrot*$kz + $zrot*$wz);

		# d
		$ncorx=$cassymdx;
		$ncory=$cassymdy;
		$ncorz=$cassymdz;
		# changement de repere
		$x=$ncorx * $ux + $ncory * $uy + $ncorz * $uz;
		$y=$ncorx * $kx + $ncory * $ky + $ncorz * $kz;
		$z=$ncorx * $wx + $ncory * $wy + $ncorz * $wz;
		# Rotation
		$xrot= $x*cos($teta) - $y*sin($teta);
		$yrot= $x*sin($teta) + $y*cos($teta);
		$zrot=$z;
		# rechangement de repere
		$cassymdx=($xrot*$ux + $yrot*$kx + $zrot*$wx);
		$cassymdy=($xrot*$uy + $yrot*$ky + $zrot*$wy);
		$cassymdz=($xrot*$uz + $yrot*$kz + $zrot*$wz);
}
else{

		# inverse les coordonnees
		$cassymax=$cassymax*-1;
		$cassymay=$cassymay*-1;
		$cassymaz=$cassymaz*-1;

		$cassymcx=$cassymcx*-1;
		$cassymcy=$cassymcy*-1;
		$cassymcz=$cassymcz*-1;

		$cassymdx=$cassymdx*-1;
		$cassymdy=$cassymdy*-1;
		$cassymdz=$cassymdz*-1;


};

## VERIFIE ANGLE

#Vecteurs
	$ux=$cassymbx-$cassymax;
	$uy=$cassymby-$cassymay;
	$uz=$cassymbz-$cassymaz;
	$nu=sqrt($ux**2+$uy**2+$uz**2);

	# vecteur v en -1 1 0
	$vx=$cassymbx+1;
	$vy=$cassymby-1;
	$vz=$cassymbz-0;
	$nv=sqrt($vx**2+$vy**2+$vz**2);

	# Angle entre les deux vecteurs a fusionner
	$cosangle=(($ux*$vx)+($uy*$vy)+($uz*$vz))/($nu * $nv);
	$testcos=1-$cosangle*$cosangle;
	$testcos=sqrt($testcos*$testcos);
	$teta=atan2(sqrt($testcos),$cosangle);
	$tetab=atan2(sqrt($testcos),$cosangle)*180.0/3.14159;
	#print"vecteur1($ux,$uy,$uz),vecteur2($vx,$vy,$vz)\n";
	#print "APRES angle(vecteur1,vecteur2) $tetab degres\n";
## FIN VERIFIE ANGLE

 	#print "$cassymax $cassymay $cassymaz\n";
 	#print "$cassymbx $cassymby $cassymbz\n";
 	#print "$cassymcx $cassymcy $cassymcz\n";
 	#print "$cassymdx $cassymdy $cassymdz\n";
	
	@tabcassymx='';
	@tabcassymy='';
	@tabcassymz='';
	
	$tabcassymx[0]=$cassymax;
	$tabcassymy[0]=$cassymay;
	$tabcassymz[0]=$cassymaz;

	$tabcassymx[1]=$cassymbx;
	$tabcassymy[1]=$cassymby;
	$tabcassymz[1]=$cassymbz;

	$tabcassymx[2]=$cassymcx;
	$tabcassymy[2]=$cassymcy;
	$tabcassymz[2]=$cassymcz;

	$tabcassymx[3]=$cassymdx;
	$tabcassymy[3]=$cassymdy;
	$tabcassymz[3]=$cassymdz;


	# cherche ou est c
	$sqtabcassymx=sqrt($tabcassymx[3]*$tabcassymx[3]);
	$sqtabcassymy=sqrt($tabcassymy[3]*$tabcassymy[3]);

	if($sqtabcassymx < $sqtabcassymy){
		$reseouz="cis" if($tabcassymx[2] < $tabcassymx[3]);
		$reseouz="trans" if($tabcassymx[2] > $tabcassymx[3]);
	}
	else{
		$reseouz="trans" if($tabcassymy[2] < $tabcassymy[3]);
		$reseouz="cis" if($tabcassymy[2] > $tabcassymy[3]);
	};

 	$reseouz;
};





