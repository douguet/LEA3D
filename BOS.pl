#!/usr/bin/perl


print"SUBROUTINE Build, Optimis and Score OK \n";

sub bos{

	local($moleculei,$nomtab)=@_;	

	unlink "mol.sdf" if(-e "mol.sdf");

	&build;#use indice $moleculei and $nomtab

	&optimis("mol.sdf") if($param{OPTIMIS} == 1 && -e "mol.sdf" && !-z "mol.sdf");
	
	if(-e "mol.sdf" && !-z "mol.sdf"){

		#chop($nbconformer = `$leaexe/NBSDF.pl mol.sdf` );
		$nbconformer=&nbsdf("mol.sdf");
		$nbconformer=~s/ //g;

		 # score (expressed in %) of conformers into @score, ranking in @ranking (in decreasing order from 1 to $nbconformers)
		 &score("mol.sdf",$nbconformer); # @ranking begins at 1

		 print "$nomtab mol[$moleculei]: best conformer is number $ranking[1] / $nbconformer with score=$score[$ranking[1]]\%\n";

	};

	#The calling subroutine can use mol.sdf, $nbconformer, @properties, @ranking and @score
};

####################################################################################
