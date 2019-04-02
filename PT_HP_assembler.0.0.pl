#!/local/bin/perl


#Ouverture des fichier d'entréé, filename = sortie de blast
#seqfile = fichier des sequences query
my $filename = $ARGV[0];
$matchfile = $ARGV[1];
#Ouverture des fichiers de sortie
$outnohit= $filename. ".nhit";
my $outfinal= $filename. ".end";
open(DAT, ">$outnohit") || die "Could not create $outnohit.\n";
open(IN, $filename) || die "Could not open $filename.\n";
open(MATCH, $matchfile) || die "Could not open $matchfile.\n";
open(END, ">$outfinal") || die "Could not create $outfinal.\n";
# utilisation des modules de perl et de bioperl, Bio::Searchio
# est le module de traitement de sortie de blast
use Bio::SearchIO;
#neconnaissance du format d'entree
    my $in = new Bio::SearchIO(-format => 'blast',
                               -file   => $filename);

#Debut du code proprement dit
#Tant que on a une query, fait la suite du programme
    while( my $result = $in->next_result ) {
#ce type de modification de variable permet de transformer les variables bioper en variables
#locales, ici c'est pour le nom de la query
    	$query_name = $result->query_name;
print "$query_name\n";
#tant que il y a des hits pour une query, continue
      while( my $hit = $result->next_hit ) {

#tant que il y a des hsp pour un hit, continue
      while( my $hsp = $hit->next_hsp ) {
#transformation des variables bioperl utilisées dans le programme en variables locales
#le nom du hit
$hit_name = $hit->name;
      print "$hit_name\n";
#le nombre d'hsp pour un hits (ne doit pas servir je crois à verifier et enlever
		 $num_hsps = $hit->num_hsps;
#le frame, attention le frame est corrige à la norme GFF,eg 1, 2, 3 resprectivement 0, 1, 2
	         $frame = $hsp->frame;
#l'orientation, tres important car par exemple -1 et 1 donne un frame de 0, l'orientation +/-
#est donnee par strand
		 $strand = $hsp->strand;
#le classement de l'hsp, on vera plus loin que seul le premier hsp est conservé
		 $rank = $hsp->rank;
#la longueur du query
		 $Qlen = $result->query_length;
#le numero d'acession seul du hit
		$accession = $hit->accession;
#la fin de l'alignement sur le query
		$Qy= $hsp->end('query');
#lancement de la premiere sous routine, se reporter plus bas à sub position_search
		$Qx= $hsp->start('query');
		$Sy= $hsp->end('subject');
		$Sx= $hsp->start('subject');
		$query_string= $hsp->query_string;
		$hit_descr= $hit->description;


seek(MATCH, 0, 0);

		while ($tp = <MATCH>){
#suprime le retour chariot

		chop($tp);
#si la ligne du fichier sequence contient signe superieur suivit du nom de la query puis fin de
#ligne, continue. attention, peut poser des probleme si on a une description dans le fichier sequence
#j'ai ete obligé de faire ca pour eviter par exemple que contig 11 soit reconnu comme contig1
 		if($tp =~ /^>/){

		if($tp =~ /\Q$hit_name/ ){

#en clair on s'est positionné au niveau de la query dans le fichier sequence
#seconde routine....
  #print END "$query_name\n";
		print "$tp\n";
		#print END "$hit_descr\n";
		&extract1_seq;
#on retient le nom de la query dans la variable $ab qui ne sera pas mise à zero dans pour le prochain
#hit
		$ab = $query_name;
		#$hit_name ="";
	 }
        }
      }
}
    }

}
& extract_nohits;



###########################################################

sub extract1_seq {

# pour rappel, on est dans la ligne juste apres le nom de la query dans le fichier sequence
#tant que la ligne ne commence pas par un signe superieur, continue
	while(($tp = <MATCH>) !~ /^>/) {
#colle la ligne dans la variable $seq, à la fin du tour de boucle (cad quand la ligne commence
#par un signe superieur, $seq sera la sequence complete
		$seqM .= $tp;
#elimine le saut de ligne
		chop($seqM);
#si fin du fichier sequence, sort de la boucle
			if(eof(MATCH))
		 {last;}
	}
&printout;
$seqM ="";
$accession= "";
last;
}

sub printout{
$Nterm = substr($seqM, 0, $Sx,);
$Cterm = substr($seqM, $Sy, );
$seqtot = join ($Cterm, $query_string, $Nterm);
$query_string =~ s/-//g;
print END ">$hit_name";
print END " $hit_descr";
print END " citrusEST ";
print END "$query_name";
print END " $Sx";
print END "-$Sy\n";
#print END "$accession\n";
print END "$Nterm";
print END "$query_string";
print END "$Cterm\n";
#print END "$seqtot\n";
}
################################

sub extract_nohits{
seek (IN, 0, 0);
print "gfgfgf";
while( my $result = $in->next_result ) {
#ce type de modification de variable permet de transformer les variables bioper en variables
#locales, ici c'est pour le nom de la query
    	my $query_name = $result->query_name;
	print "$query_name\n";
#si il n'y a pas de hit pour une query, alors imprime dans le fichier .DAT le nom de la query
	if($result->next_hit eq ""){
	print DAT "$query_name\n";
}

}

}





