#!/usr/bin/perl


print"SUBROUTINE OPTIMIS OK \n";

sub optimis{

	local($file)=@_;	
	
	$fileco=$file;
	$fileco=~s/\.sdf$/co\.sdf/;
	unlink $fileco;
	#system("touch $fileco"); #must be initialized !!!

	$platform=$^O;

	$param{NBCONF}=1 if($param{NBCONF} eq "" || $param{NBCONF} < 1);
	
	#$convertorfrog=$ENV{LEA_FROG};
	#$convertorcorina=$ENV{LEA_CORINA};
	#$convertorrdkit=$ENV{LEA_RDKIT};
	#$convertorrdkitmin=$ENV{LEA_RDKITmin};
	$convertorfrog="";
	$convertorcorina="";
	$convertorrdkitmin="";
	$convertorrdkit=$convertor;

	if($convertor eq $convertorcorina){

		#******************************************************
		#	CORINA: input= $file output= $fileco (multisdf)
		#*******************************************************

		unlink "corina.trc" if(-e "corina.trc");

		if($param{NBCONF}==1){
			system("$convertor -i t=sdf,sdfi2n=ID -o t=sdf -t n -d wh $param{WORKDIR}/$file $param{WORKDIR}/$fileco");
		}
		else{
			system("$convertor -i t=sdf,sdfi2n=ID -o t=sdf -t n -d wh,rc,mc=$param{NBCONF},r2d,flapn,sc $param{WORKDIR}/$file $param{WORKDIR}/$fileco");
		};

		print "WARNING : SEE corina.trc\n" if(-e "corina.trc" && ! -z "corina.trc");

		#*******************************************************	
	}
	elsif($convertor eq $convertorrdkit){
		#**************************************************
		#	RDKIT
		#**************************************************
		if($platform eq "linux" || $platform eq "darwin"){
			system("$convertor $param{WORKDIR}/$file $param{NBCONF} $param{WORKDIR}/$fileco");
		}
		else{#windows
			system("python $convertor $param{WORKDIR}/$file $param{NBCONF} $param{WORKDIR}/$fileco");
		};
		#system(" $convertor $param{WORKDIR}/$file $param{NBCONF} $param{WORKDIR}/$fileco");
	}
	elsif($convertor eq $convertorrdkitmin){
		#**************************************************
		#       RDKITmin
		#**************************************************
		# 1 minimized conformer only
		system("$convertor $param{WORKDIR}/$file $param{WORKDIR}/$fileco");
	}	
	else{
		#**************************************************
		#	FROG
		#**************************************************
		# use the program frog in $leaexe
		#frog output is <>3D.sdf

		$frog3d=$file;
		$frog3d=~s/\.sdf$/3D\.sdf/;
		unlink "$frog3d";
		if($param{NBCONF}==1){
			system("$leaexe/frog $param{WORKDIR}/$file sdf 1");
		}
		else{
			system("$leaexe/frog $param{WORKDIR}/$file sdf $param{NBCONF}");
		};
		if(-e "$frog3d" && !-z "$frog3d"){
			rename "$frog3d", "$fileco";
		};
	};	


	if(! -e $fileco || -z $fileco){
		print "WARNING: optimization failed, use $file \n";
		$memo_remark[$memoi]=$memo_remark[$memoi]."optimization failed (use the non optimized sdf file)\n" if((!-e $fileco || -z $fileco) && $memo_remark[$memoi] !~/optimization failed/);
	}
	else{
		rename $fileco, $file;
	};	
	
	#chop($cms1 = `$leaexe/NBSDF.pl $file` );
	$cms1=&nbsdf($file);
	print "$file : $cms1 molecules\n" if($param{VERBOSITY} >= 1);

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
