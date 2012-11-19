#!/usr/bin/perl 
use strict;
use Getopt::Std;
use List::Util qw ( min max sum );

use lib '/home/tomfy/MHarrisonProject/lib';
use Overlap;
use Orthologger;

use lib '/home/tomfy/cxgn/cxgn-corelibs/lib';
use CXGN::Phylo::File;
use CXGN::Phylo::Parser;

use Devel::Cycle;
# find_cycle($test); # to find circular refs in $test
         

$ENV{PERL_DEBUG_MSTATS} = 2;

my $N = shift || 100;
my $gene_tree_newick = 
'(((((Bradi1g72517.1[species=Brachypodium_distachyon]:0.084953,((Sb06g001440.1[
species=Sorghum_bicolor]:0.053689,GRMZM2G167658_P01[species=Zea_mays]:0.044051)[
speciation=1]:0.049117,LOC_Os03g08380.1[species=Oryza_sativa]:0.093295)[speciati
on=1]:0.006628)[speciation=0]:0.140355,(((Glyma16g01350.1[species=Glycine_max]:0
.101234,IMGA_Medtr8g022270.1[species=Medicago_truncatula]:0.084722)[speciation=1
]:0.051900,((evm.model.supercontig_25.158[species=Carica_papaya]:0.561998,POPTR_
0001s16560.1[species=Populus_trichocarpa]:0.011892)[speciation=1]:0.050806,X_297
27.m000475[species=Ricinus_communis]:0.054235)[speciation=0]:0.050869)[speciatio
n=0]:0.019812,GSVIVT01016706001[species=Vitis_vinifera]:0.242449)[speciation=1]:
0.078188)[speciation=1]:0.361309,(((((LOC_Os08g45030.1[species=Oryza_sativa]:0.0
73065,(Sb07g023730.1[species=Sorghum_bicolor]:0.011580,GRMZM2G315375_P01[species
=Zea_mays]:0.018349)[speciation=1]:0.048954)[speciation=1]:0.011967,Bradi3g12627
.1[species=Brachypodium_distachyon]:0.082723)[speciation=0]:0.056375,((AT2G36910.1[species=Arabidopsis_thaliana]:0.045687,(((POPTR_0006s12590.1[species=Populus_tr
ichocarpa]:0.025476,POPTR_0016s09680.1[species=Populus_trichocarpa]:0.015892)[speciation=0]:0.012899,((((GSVIVT01033645001[species=Vitis_vinifera]:0.057919,Glyma1
3g20530.1[species=Glycine_max]:0.223115)[speciation=1]:0.351768,Glyma10g06220.1[species=Glycine_max]:0.0001)[speciation=0]:0.053756,Glyma19g36820.1[species=Glycin
e_max]:0.0001)[speciation=0]:0.012752,Glyma03g34080.1[species=Glycine_max]:0.009798)[speciation=0]:0.044684)[speciation=0]:0.001994,X_30078.m002286[species=Ricinu
s_communis]:0.033727)[speciation=0]:0.016653)[speciation=0]:0.013972,Solyc09g008240.2.1[species=Solanum_lycopersicum]:0.063860)[speciation=1]:0.054056)[speciation
=1]:0.247162,((((Bradi5g12307.1[species=Brachypodium_distachyon]:0.021628,(LOC_Os04g38570.1[species=Oryza_sativa]:0.014406,(GRMZM2G125424_P04[species=Zea_mays]:0.
008203,Sb06g018860.1[species=Sorghum_bicolor]:0.004471)[speciation=1]:0.017528)[speciation=1]:0.005279)[speciation=0]:0.050064,((((LOC_Os04g38570.2[species=Oryza_
sativa]:0.082651,GRMZM2G072850_P02[species=Zea_mays]:0.081078)[speciation=1]:0.145575,GRMZM2G085236_P02[species=Zea_mays]:0.019828)[speciation=0]:0.063454,((GRMZM
2G072850_P01[species=Zea_mays]:0.008045,GRMZM2G085236_P01[species=Zea_mays]:0.008905)[speciation=0]:0.001653,Sb06g030350.1[species=Sorghum_bicolor]:0.005113)[spec
iation=1]:0.003957)[speciation=0]:0.017051,Bradi5g23600.1[species=Brachypodium_distachyon]:0.022423)[speciation=0]:0.044032)[speciation=0]:0.019204,(AT3G28860.1[s
pecies=Arabidopsis_thaliana]:0.056394,(Solyc02g087870.2.1[species=Solanum_lycopersicum]:0.036583,((evm.model.supercontig_14.14[species=Carica_papaya]:0.030364,((X
_29822.m003496[species=Ricinus_communis]:0.017461,POPTR_0017s11750.1[species=Populus_trichocarpa]:0.018605)[speciation=1]:0.007765,(Glyma13g05300.1[species=Glycin
e_max]:0.001317,Glyma19g02520.1[species=Glycine_max]:0.002883)[speciation=0]:0.029614)[speciation=1]:0.002277)[speciation=1]:0.005123,GSVIVT01011381001[species=Vi
tis_vinifera]:0.154677)[speciation=1]:0.004627)[speciation=1]:0.012197)[speciation=0]:0.035723)[speciation=1]:0.255196,((((Bradi3g52220.1[species=Brachypodium_dis
tachyon]:0.051750,(LOC_Os02g46680.1[species=Oryza_sativa]:0.049241,(GRMZM2G004748_P01[species=Zea_mays]:0.027222,Sb04g031170.1[species=Sorghum_bicolor]:0.019569)[
speciation=1]:0.016525)[speciation=1]:0.006971)[speciation=0]:0.095426,((AT1G10680.1[species=Arabidopsis_thaliana]:0.117913,AT4G25960.1[species=Arabidopsis_thalia
na]:0.085843)[speciation=0]:0.037122,(Solyc08g076720.2.1[species=Solanum_lycopersicum]:0.108681,((POPTR_0001s02220.1[species=Populus_trichocarpa]:0.111959,(((X_30
054.m000810[species=Ricinus_communis]:0.013004,POPTR_0003s09320.2[species=Populus_trichocarpa]:0.438621)[speciation=1]:0.111967,POPTR_0003s09320.1[species=Populus
_trichocarpa]:0.0001)[speciation=0]:0.053813,POPTR_0001s02200.1[species=Populus_trichocarpa]:0.040292)[speciation=0]:0.054962)[speciation=0]:0.026803,((evm.model.
supercontig_90.67[species=Carica_papaya]:0.072710,GSVIVT01013125001[species=Vitis_vinifera]:0.097435)[speciation=1]:0.009214,((IMGA_Medtr5g029750.1[species=Medica
go_truncatula]:0.044084,(Glyma01g02060.1[species=Glycine_max]:0.006103,Glyma09g33880.1[species=Glycine_max]:0.004887)[speciation=0]:0.043860)[speciation=1]:0.0216
15,Glyma08g36450.1[species=Glycine_max]:0.264192)[speciation=0]:0.023583)[speciation=0]:0.004439)[speciation=0]:0.012825)[speciation=1]:0.014366)[speciation=0]:0.
020906)[speciation=1]:0.210399,(((AT1G28010.1[species=Arabidopsis_thaliana]:0.067506,AT1G27940.1[species=Arabidopsis_thaliana]:0.061446)[speciation=0]:0.133524,((
POPTR_0002s02110.1[species=Populus_trichocarpa]:0.138425,(X_30170.m014384[species=Ricinus_communis]:0.099901,GSVIVT01009946001[species=Vitis_vinifera]:0.403313)[s
peciation=1]:0.013873)[speciation=0]:0.029678,(IMGA_Medtr1g025560.1[species=Medicago_truncatula]:0.155051,(Glyma14g40280.1[species=Glycine_max]:0.068215,Glyma17g3
7860.1[species=Glycine_max]:0.025615)[speciation=0]:0.065247)[speciation=1]:0.055277)[speciation=0]:0.040112)[speciation=0]:0.188653,((jgi_Selmo1_230134[species=S
elaginella]:0.031172,jgi_Selmo1_236608[species=Selaginella]:0.050443)[speciation=0]:0.309863,((jgi_Selmo1_105467[species=Selaginella]:0.008053,jgi_Selmo1_148837[s
pecies=Selaginella]:0.014924)[speciation=0]:0.178585,(jgi_Selmo1_230688[species=Selaginella]:0.054324,jgi_Selmo1_117529[species=Selaginella]:0.0001)[speciation=0]
:0.187892)[speciation=0]:0.085036)[speciation=0]:0.012789)[speciation=1]:0.019635)[speciation=0]:0.022028,((jgi_Selmo1_410515[species=Selaginella]:0.007408,jgi_Se
lmo1_138662[species=Selaginella]:0.039342)[speciation=0]:0.357277,(jgi_Selmo1_92485[species=Selaginella]:0.012331,jgi_Selmo1_114581[species=Selaginella]:0.013289)
[speciation=0]:0.286680)[speciation=0]:0.035912)[speciation=0]:0.004513)[speciation=0]:0.0001)[speciation=0]:0.049965,((((jgi_Selmo1_84555[species=Selaginella]:0.
005219,jgi_Selmo1_129540[species=Selaginella]:0.003306)[speciation=0]:0.464418,(jgi_Selmo1_86998[species=Selaginella]:0.015196,jgi_Selmo1_75690[species=Selaginell
a]:0.013109)[speciation=0]:0.494939)[speciation=0]:0.015547,((((LOC_Os08g05690.1[species=Oryza_sativa]:0.062840,Bradi3g17020.1[species=Brachypodium_distachyon]:0.
071615)[speciation=1]:0.024601,(Sb07g003510.1[species=Sorghum_bicolor]:0.033130,GRMZM2G333183_P01[species=Zea_mays]:0.034693)[speciation=1]:0.042529)[speciation=1
]:0.111295,(Bradi3g17010.1[species=Brachypodium_distachyon]:0.114720,((GRMZM2G441722_P01[species=Zea_mays]:0.169627,Sb07g003520.1[species=Sorghum_bicolor]:0.0001)
[speciation=1]:0.102827,LOC_Os08g05710.1[species=Oryza_sativa]:0.108029)[speciation=1]:0.021166)[speciation=0]:0.120020)[speciation=0]:0.065184,(((((POPTR_0015s00
250.1[species=Populus_trichocarpa]:0.133082,X_30098.m001722[species=Ricinus_communis]:0.126658)[speciation=1]:0.009542,GSVIVT01038687001[species=Vitis_vinifera]:0
.225975)[speciation=1]:0.021811,evm.model.supercontig_3.192[species=Carica_papaya]:0.204476)[speciation=0]:0.017473,(IMGA_Medtr3g093430.1[species=Medicago_truncat
ula]:0.076671,Glyma06g14450.1[species=Glycine_max]:0.079281)[speciation=1]:0.086749)[speciation=0]:0.022673,Solyc06g072960.1.1[species=Solanum_lycopersicum]:0.374
125)[speciation=1]:0.071074)[speciation=1]:0.352057)[speciation=1]:0.017819,((((((((((GSVIVT01025040001[species=Vitis_vinifera]:0.376890,GSVIVT01016617001[species
=Vitis_vinifera]:0.010049)[speciation=0]:0.074386,X_28180.m000390[species=Ricinus_communis]:0.002556)[speciation=1]:0.105207,((evm.model.supercontig_1.347[species
=Carica_papaya]:0.603456,Glyma10g43700.1[species=Glycine_max]:0.0001)[speciation=1]:0.042849,Glyma20g38380.1[species=Glycine_max]:0.0001)[speciation=0]:0.033345)[
speciation=0]:0.016208,POPTR_0010s21720.1[species=Populus_trichocarpa]:0.011181)[speciation=0]:0.008968,POPTR_0008s05020.1[species=Populus_trichocarpa]:0.016802)[
speciation=0]:0.028848,Solyc04g010310.2.1[species=Solanum_lycopersicum]:0.087783)[speciation=1]:0.010457,(Glyma02g10530.1[species=Glycine_max]:0.008474,Glyma18g52
350.1[species=Glycine_max]:0.011966)[speciation=0]:0.048697)[speciation=0]:0.006728,(AT2G39480.1[species=Arabidopsis_thaliana]:0.044322,AT3G55320.1[species=Arabid
opsis_thaliana]:0.033617)[speciation=0]:0.036223)[speciation=0]:0.040076,(((Bradi1g66290.1[species=Brachypodium_distachyon]:0.060919,(GRMZM5G891159_P01[species=Ze
a_mays]:0.011856,Sb01g039110.1[species=Sorghum_bicolor]:0.010314)[speciation=1]:0.013923)[speciation=1]:0.009766,LOC_Os03g17180.1[species=Oryza_sativa]:0.015982)[
speciation=0]:0.057236,((Bradi2g62540.1[species=Brachypodium_distachyon]:0.027052,LOC_Os01g74470.1[species=Oryza_sativa]:0.014423)[speciation=1]:0.006428,Sb03g047
490.1[species=Sorghum_bicolor]:0.021011)[speciation=1]:0.063109)[speciation=0]:0.006621)[speciation=1]:0.139668,(jgi_Selmo1_169170[species=Selaginella]:0.000119,j
gi_Selmo1_139408[species=Selaginella]:0.000721)[speciation=0]:0.248346)[speciation=1]:0.361856)[speciation=0]:0.008863)[speciation=0]:0.038646)[speciation=0]:0.00
01,((((((GSVIVT01007586001[species=Vitis_vinifera]:0.353203,(POPTR_0006s07370.1[species=Populus_trichocarpa]:0.121459,(POPTR_0018s13810.1[species=Populus_trichoca
rpa]:0.062299,(POPTR_0006s07390.1[species=Populus_trichocarpa]:0.143030,POPTR_0018s13820.1[species=Populus_trichocarpa]:0.636231)[speciation=0]:0.086438)[speciati
on=0]:0.069687)[speciation=0]:0.085187)[speciation=1]:0.157604,(((GSVIVT01014625001[species=Vitis_vinifera]:0.302629,GSVIVT01000580001[species=Vitis_vinifera]:0.2
88478)[speciation=0]:0.129244,GSVIVT01036801001[species=Vitis_vinifera]:0.685686)[speciation=0]:0.177576,(Solyc07g018130.1.1[species=Solanum_lycopersicum]:0.26501
5,((Solyc07g064120.1.1[species=Solanum_lycopersicum]:0.188940,((IMGA_Medtr4g081190.1[species=Medicago_truncatula]:0.088590,Glyma06g42040.1[species=Glycine_max]:0.
141638)[speciation=1]:0.073485,(evm.model.supercontig_97.65[species=Carica_papaya]:0.134442,((X_28166.m001066[species=Ricinus_communis]:0.117475,POPTR_0001s44320.
1[species=Populus_trichocarpa]:0.096763)[speciation=1]:0.022969,POPTR_0011s13720.1[species=Populus_trichocarpa]:0.129953)[speciation=0]:0.027626)[speciation=1]:0.
024976)[speciation=0]:0.032017)[speciation=1]:0.062600,(IMGA_Medtr3g086430.1[species=Medicago_truncatula]:0.275999,X_29929.m004786[species=Ricinus_communis]:0.206
214)[speciation=1]:0.022960)[speciation=0]:0.018182)[speciation=0]:0.037662)[speciation=0]:0.021262)[speciation=0]:0.030865,((((IMGA_Medtr7g051100.1[species=Medic
ago_truncatula]:0.155134,Glyma08g45660.1[species=Glycine_max]:0.154058)[speciation=1]:0.045876,((((IMGA_Medtr6g008820.1[species=Medicago_truncatula]:0.139154,IMGA
_Medtr6g008800.1[species=Medicago_truncatula]:0.115105)[speciation=0]:0.012036,(IMGA_Medtr6g009110.1[species=Medicago_truncatula]:0.081875,(IMGA_Medtr6g009030.1[s
pecies=Medicago_truncatula]:0.061480,(IMGA_Medtr6g009150.1[species=Medicago_truncatula]:0.060312,IMGA_Medtr6g009200.1[species=Medicago_truncatula]:0.043982)[speci
ation=0]:0.016333)[speciation=0]:0.011367)[speciation=0]:0.009973)[speciation=0]:0.020512,Glyma19g01940.1[species=Glycine_max]:0.109316)[speciation=1]:0.032996,(G
lyma19g01970.1[species=Glycine_max]:0.098720,Glyma19g01980.1[species=Glycine_max]:0.102254)[speciation=0]:0.102480)[speciation=0]:0.020487)[speciation=0]:0.019856
,((((POPTR_0017s12120.1[species=Populus_trichocarpa]:0.150463,(POPTR_0017s11030.1[species=Populus_trichocarpa]:0.101486,GSVIVT01032578001[species=Vitis_vinifera]:
0.181435)[speciation=1]:0.011850)[speciation=0]:0.012595,evm.model.supercontig_14.71[species=Carica_papaya]:0.109469)[speciation=0]:0.018998,Solyc02g087410.2.1[sp
ecies=Solanum_lycopersicum]:0.135591)[speciation=1]:0.023821,(((AT3G28415.1[species=Arabidopsis_thaliana]:0.134413,AT3G28380.1[species=Arabidopsis_thaliana]:0.111
678)[speciation=0]:0.002000,(AT3G28360.1[species=Arabidopsis_thaliana]:0.121651,AT3G28390.1[species=Arabidopsis_thaliana]:0.094819)[speciation=0]:0.004680)[specia
tion=0]:0.045149,AT3G28345.1[species=Arabidopsis_thaliana]:0.113801)[speciation=0]:0.051393)[speciation=0]:0.017332)[speciation=0]:0.049198,((Sb04g006090.1[specie
s=Sorghum_bicolor]:0.165375,(((Sb04g006100.1[species=Sorghum_bicolor]:0.002999,Sb04g006087.1[species=Sorghum_bicolor]:0.499046)[speciation=0]:0.030077,LOC_Os02g09
720.1[species=Oryza_sativa]:0.027853)[speciation=1]:0.018752,Bradi3g06577.1[species=Brachypodium_distachyon]:0.049909)[speciation=0]:0.039954)[speciation=0]:0.028
029,(GRMZM5G843192_P01[species=Zea_mays]:0.039986,Sb04g022480.1[species=Sorghum_bicolor]:0.059495)[speciation=1]:0.114473)[speciation=0]:0.097651)[speciation=1]:0
.070074)[speciation=0]:0.018730,((GRMZM2G111462_P01[species=Zea_mays]:0.057560,Sb06g020350.1[species=Sorghum_bicolor]:0.070276)[speciation=1]:0.052215,(Bradi5g136
17.1[species=Brachypodium_distachyon]:0.112478,LOC_Os04g40570.1[species=Oryza_sativa]:0.117497)[speciation=1]:0.012983)[speciation=1]:0.288107)[speciation=0]:0.03
8877,(((Bradi2g48610.1[species=Brachypodium_distachyon]:0.052401,LOC_Os01g52550.1[species=Oryza_sativa]:0.048494)[speciation=1]:0.007695,(GRMZM2G025860_P01[specie
s=Zea_mays]:0.016350,Sb03g033290.1[species=Sorghum_bicolor]:0.011081)[speciation=1]:0.033018)[speciation=1]:0.151294,(Solyc03g093650.2.1[species=Solanum_lycopersi
cum]:0.153219,((Glyma01g01160.1[species=Glycine_max]:0.046810,Glyma16g08480.1[species=Glycine_max]:0.045116)[speciation=0]:0.087490,(evm.model.supercontig_40.105[
species=Carica_papaya]:0.115705,(GSVIVT01015306001[species=Vitis_vinifera]:0.297945,(X_29912.m005427[species=Ricinus_communis]:0.065890,POPTR_0018s09420.1[species
=Populus_trichocarpa]:0.070587)[speciation=1]:0.008718)[speciation=1]:0.013718)[speciation=0]:0.023936)[speciation=0]:0.032830)[speciation=1]:0.046321)[speciation
=1]:0.214585)[speciation=0]:0.056213,(((jgi_Selmo1_229984[species=Selaginella]:0.051669,jgi_Selmo1_117838[species=Selaginella]:0.018215)[speciation=0]:0.265240,((
jgi_Selmo1_99885[species=Selaginella]:0.003178,jgi_Selmo1_103646[species=Selaginella]:0.006190)[speciation=0]:0.247511,((jgi_Selmo1_235518[species=Selaginella]:0.
0001,jgi_Selmo1_419451[species=Selaginella]:0.521301)[speciation=0]:0.109197,jgi_Selmo1_233433[species=Selaginella]:0.0001)[speciation=0]:0.363613)[speciation=0]:
0.080917)[speciation=0]:0.055655,((jgi_Selmo1_429505[species=Selaginella]:0.005919,jgi_Selmo1_231714[species=Selaginella]:0.089720)[speciation=0]:0.567648,(jgi_Se
lmo1_408755[species=Selaginella]:0.000249,jgi_Selmo1_421121[species=Selaginella]:0.005698)[speciation=0]:0.432679)[speciation=0]:0.057425)[speciation=0]:0.033576)
[speciation=1]:0.071346)[speciation=0]:0.081453,((((jgi_Selmo1_169182[species=Selaginella]:0.000204,jgi_Selmo1_154740[species=Selaginella]:0.000636)[speciation=0]
:0.227284,(jgi_Selmo1_176522[species=Selaginella]:0.010629,jgi_Selmo1_177681[species=Selaginella]:0.011484)[speciation=0]:0.263287)[speciation=0]:0.006826,((jgi_S
elmo1_74892[species=Selaginella]:0.004871,jgi_Selmo1_87743[species=Selaginella]:0.006952)[speciation=0]:0.268400,((jgi_Selmo1_123936[species=Selaginella]:0.010919
,jgi_Selmo1_184091[species=Selaginella]:0.010331)[speciation=0]:0.200066,(jgi_Selmo1_123915[species=Selaginella]:0.001401,jgi_Selmo1_184079[species=Selaginella]:0
.003642)[speciation=0]:0.200771)[speciation=0]:0.033550)[speciation=0]:0.007021)[speciation=0]:0.061304,(((Solyc11g067310.1.1[species=Solanum_lycopersicum]:0.1445
86,Solyc11g067300.1.1[species=Solanum_lycopersicum]:0.136191)[speciation=0]:0.217423,(((Bradi2g36897.1[species=Brachypodium_distachyon]:0.069238,Sb09g002940.1[spe
cies=Sorghum_bicolor]:0.074710)[speciation=1]:0.019980,(LOC_Os01g18670.1[species=Oryza_sativa]:0.040490,(((Sb03g011860.1[species=Sorghum_bicolor]:0.009822,GRMZM2G
085111_P02[species=Zea_mays]:0.016623)[speciation=1]:0.004897,(GRMZM2G119894_P01[species=Zea_mays]:0.0001,GRMZM2G119894_P03[species=Zea_mays]:0.0001)[speciation=0
]:0.019285)[speciation=0]:0.020887,Bradi2g11210.2[species=Brachypodium_distachyon]:0.035609)[speciation=1]:0.003395)[speciation=0]:0.032610)[speciation=0]:0.13397
6,((((((((X_30076.m004556[species=Ricinus_communis]:0.117832,(POPTR_0002s18850.1[species=Populus_trichocarpa]:0.049919,POPTR_0014s10860.1[species=Populus_trichoca
rpa]:0.044781)[speciation=0]:0.061420)[speciation=1]:0.026854,((((Glyma02g01100.1[species=Glycine_max]:0.014415,Glyma10g27790.1[species=Glycine_max]:0.020771)[spe
ciation=0]:0.024457,IMGA_Medtr1g086080.1[species=Medicago_truncatula]:0.040357)[speciation=1]:0.028875,Glyma03g38300.1[species=Glycine_max]:0.096374)[speciation=0
]:0.040767,(((X_30076.m004557[species=Ricinus_communis]:0.029069,X_30076.m004558[species=Ricinus_communis]:0.027535)[speciation=0]:0.069497,((POPTR_0014s10870.1[s
pecies=Populus_trichocarpa]:0.099113,((POPTR_0002s18860.2[species=Populus_trichocarpa]:0.003277,POPTR_0002s18860.1[species=Populus_trichocarpa]:0.000115)[speciati
on=0]:0.053949,(POPTR_0014s10880.1[species=Populus_trichocarpa]:0.0001,POPTR_0014s10880.2[species=Populus_trichocarpa]:0.0001)[speciation=0]:0.041922)[speciation=
0]:0.027478)[speciation=0]:0.007279,X_30076.m004559[species=Ricinus_communis]:0.125074)[speciation=1]:0.009798)[speciation=0]:0.008396,(GSVIVT01028256001[species=
Vitis_vinifera]:0.137475,evm.model.supercontig_57.7[species=Carica_papaya]:0.208226)[speciation=1]:0.070441)[speciation=0]:0.011345)[speciation=0]:0.005137)[speci
ation=0]:0.014470,(Solyc06g009290.2.1[species=Solanum_lycopersicum]:0.049509,Solyc06g009280.1.1[species=Solanum_lycopersicum]:0.070398)[speciation=0]:0.088317)[sp
eciation=1]:0.007771,(((AT1G02520.1[species=Arabidopsis_thaliana]:0.025530,AT1G02530.1[species=Arabidopsis_thaliana]:0.051271)[speciation=0]:0.062515,(AT4G01820.1
[species=Arabidopsis_thaliana]:0.077598,AT4G01830.1[species=Arabidopsis_thaliana]:0.072482)[speciation=0]:0.058656)[speciation=0]:0.059223,(AT2G47000.1[species=Ar
abidopsis_thaliana]:0.083346,AT3G62150.1[species=Arabidopsis_thaliana]:0.085427)[speciation=0]:0.062974)[speciation=0]:0.011966)[speciation=0]:0.007041,((((IMGA_M
edtr2g018350.1[species=Medicago_truncatula]:0.319100,IMGA_Medtr1g086150.1[species=Medicago_truncatula]:0.306602)[speciation=0]:0.073331,(((((((GRMZM2G014089_P02[s
pecies=Zea_mays]:0.485895,((evm.model.supercontig_102.36[species=Carica_papaya]:0.113011,GRMZM2G315375_P02[species=Zea_mays]:0.181973)[speciation=1]:0.218182,IMGA
_Medtr8g066690.1[species=Medicago_truncatula]:0.472731)[speciation=0]:0.158861)[speciation=0]:0.239516,((X_29822.m003426[species=Ricinus_communis]:0.229714,IMGA_M
edtr6g009060.1[species=Medicago_truncatula]:0.258620)[speciation=1]:0.318995,IMGA_Medtr6g009090.1[species=Medicago_truncatula]:0.0001)[speciation=0]:0.265760)[spe
ciation=0]:0.277988,IMGA_Medtr7g102070.1[species=Medicago_truncatula]:0.253126)[speciation=0]:0.118974,LOC_Os04g54930.1[species=Oryza_sativa]:0.293475)[speciation
=0]:0.061224,(GRMZM2G401769_P01[species=Zea_mays]:0.031417,GRMZM2G049351_P01[species=Zea_mays]:0.463266)[speciation=0]:0.535028)[speciation=0]:0.133522,(GRMZM5G83
3207_P01[species=Zea_mays]:1.148219,GRMZM5G843537_P01[species=Zea_mays]:0.078264)[speciation=0]:0.317539)[speciation=0]:0.205285,Glyma18g24280.1[species=Glycine_m
ax]:0.492519)[speciation=0]:0.246825)[speciation=0]:0.238230,evm.model.supercontig_57.8[species=Carica_papaya]:0.271540)[speciation=0]:0.262092,(POPTR_0010s00540.
1[species=Populus_trichocarpa]:0.144421,X_29889.m003409[species=Ricinus_communis]:0.200409)[speciation=1]:0.031873)[speciation=0]:0.015313)[speciation=0]:0.015159
,(GSVIVT01017696001[species=Vitis_vinifera]:0.237976,X_29889.m003408[species=Ricinus_communis]:0.211765)[speciation=1]:0.013562)[speciation=0]:0.002592,((((((Glym
a13g17910.1[species=Glycine_max]:0.051749,(Glyma17g04600.1[species=Glycine_max]:0.237407,Glyma13g17920.1[species=Glycine_max]:0.020378)[speciation=0]:0.036221)[sp
eciation=0]:0.034668,(Glyma17g04590.1[species=Glycine_max]:0.021597,(Glyma13g17930.1[species=Glycine_max]:0.0001,Glyma13g17930.2[species=Glycine_max]:0.136922)[sp
eciation=0]:0.064442)[speciation=0]:0.044045)[speciation=0]:0.005456,((IMGA_Medtr4g077930.1[species=Medicago_truncatula]:0.0001,IMGA_Medtr4g078030.1[species=Medic
ago_truncatula]:0.0001)[speciation=0]:0.066537,(IMGA_Medtr4g123990.1[species=Medicago_truncatula]:0.031792,IMGA_Medtr4g124000.1[species=Medicago_truncatula]:0.036
680)[speciation=0]:0.023091)[speciation=0]:0.053769)[speciation=1]:0.001381,IMGA_Medtr4g124040.1[species=Medicago_truncatula]:0.065313)[speciation=0]:0.057048,(((
Glyma17g04610.1[species=Glycine_max]:0.043021,Glyma13g17890.1[species=Glycine_max]:0.137404)[speciation=0]:0.039635,(Glyma13g17880.1[species=Glycine_max]:0.444462
,Glyma17g04620.1[species=Glycine_max]:0.055997)[speciation=0]:0.075224)[speciation=0]:0.020590,((IMGA_Medtr4g124050.1[species=Medicago_truncatula]:0.0001,IMGA_Med
tr4g124050.2[species=Medicago_truncatula]:0.141829)[speciation=0]:0.096395,IMGA_Medtr3g107800.1[species=Medicago_truncatula]:0.070417)[speciation=0]:0.039688)[spe
ciation=1]:0.080524)[speciation=0]:0.005530,IMGA_Medtr3g080220.1[species=Medicago_truncatula]:0.182598)[speciation=0]:0.032881)[speciation=0]:0.004886,(Solyc12g09
8840.1.1[species=Solanum_lycopersicum]:0.051913,Solyc12g098870.1.1[species=Solanum_lycopersicum]:0.057226)[speciation=0]:0.214317)[speciation=0]:0.013300)[speciat
ion=0]:0.017567)[speciation=0]:0.043302,((X_29581.m000253[species=Ricinus_communis]:0.343455,(Solyc03g005860.2.1[species=Solanum_lycopersicum]:0.260705,((((X_2969
3.m002086[species=Ricinus_communis]:0.150527,(POPTR_0001s34280.1[species=Populus_trichocarpa]:0.117888,X_29693.m002084[species=Ricinus_communis]:0.122419)[speciat
ion=1]:0.005591)[speciation=0]:0.033380,((((IMGA_Medtr2g018320.1[species=Medicago_truncatula]:0.092681,Glyma15g09680.1[species=Glycine_max]:0.254779)[speciation=1
]:0.012754,(IMGA_Medtr6g078080.1[species=Medicago_truncatula]:0.083481,IMGA_Medtr2g018530.1[species=Medicago_truncatula]:0.068653)[speciation=0]:0.050044)[speciat
ion=0]:0.013403,Glyma13g29380.1[species=Glycine_max]:0.129636)[speciation=0]:0.041387,(GSVIVT01021366001[species=Vitis_vinifera]:0.185189,GSVIVT01021365001[specie
s=Vitis_vinifera]:0.187370)[speciation=0]:0.020124)[speciation=1]:0.007069)[speciation=0]:0.003243,(AT5G46540.1[species=Arabidopsis_thaliana]:0.176644,AT4G18050.1
[species=Arabidopsis_thaliana]:0.133487)[speciation=0]:0.053481)[speciation=0]:0.002310,(Solyc02g071350.2.1[species=Solanum_lycopersicum]:0.084705,Solyc02g071340.
1.1[species=Solanum_lycopersicum]:0.122872)[speciation=0]:0.091174)[speciation=1]:0.025263)[speciation=0]:0.014111)[speciation=0]:0.034410,((((((((LOC_Os01g35030.
1[species=Oryza_sativa]:0.085439,((LOC_Os01g34970.1[species=Oryza_sativa]:0.0001,LOC_Os01g34970.2[species=Oryza_sativa]:0.316556)[speciation=0]:0.624009,Bradi2g20
177.1[species=Brachypodium_distachyon]:0.0001)[speciation=1]:0.127005)[speciation=0]:0.056050,Sb03g023740.1[species=Sorghum_bicolor]:0.093541)[speciation=1]:0.038
607,GRMZM2G086730_P01[species=Zea_mays]:0.058223)[speciation=0]:0.023696,Sb02g019540.1[species=Sorghum_bicolor]:0.098601)[speciation=0]:0.016819,Bradi2g20170.1[sp
ecies=Brachypodium_distachyon]:0.132948)[speciation=0]:0.135060,(LOC_Os05g47490.1[species=Oryza_sativa]:0.108365,(GRMZM2G082385_P01[species=Zea_mays]:0.054331,Sb0
9g027320.1[species=Sorghum_bicolor]:0.062624)[speciation=1]:0.065491)[speciation=1]:0.121514)[speciation=0]:0.016009,(((((Bradi3g20045.1[species=Brachypodium_dist
achyon]:0.143487,(Bradi2g47410.1[species=Brachypodium_distachyon]:0.058237,Bradi2g47400.1[species=Brachypodium_distachyon]:0.037417)[speciation=0]:0.027279)[speci
ation=0]:0.020582,(LOC_Os01g50160.1[species=Oryza_sativa]:0.063734,Sb03g032030.1[species=Sorghum_bicolor]:0.081232)[speciation=1]:0.023959)[speciation=0]:0.058879
,(((Bradi2g47337.1[species=Brachypodium_distachyon]:0.054271,(Bradi2g47347.1[species=Brachypodium_distachyon]:0.085568,Bradi2g47360.1[species=Brachypodium_distach
yon]:0.033354)[speciation=0]:0.014216)[speciation=0]:0.025263,(GRMZM2G014089_P01[species=Zea_mays]:0.345093,Sb03g032000.1[species=Sorghum_bicolor]:0.023722)[speci
ation=1]:0.097993)[speciation=1]:0.023098,LOC_Os01g50100.1[species=Oryza_sativa]:0.142659)[speciation=0]:0.087642)[speciation=0]:0.012702,(LOC_Os01g50080.1[specie
s=Oryza_sativa]:0.180214,((AC233882.1_FGP003[species=Zea_mays]:0.061656,Sb03g031990.1[species=Sorghum_bicolor]:0.045542)[speciation=1]:0.064784,Bradi2g47330.1[spe
cies=Brachypodium_distachyon]:0.107649)[speciation=1]:0.003684)[speciation=0]:0.112582)[speciation=0]:0.017452,(((Bradi2g17720.1[species=Brachypodium_distachyon]:
0.076265,Bradi2g17710.1[species=Brachypodium_distachyon]:0.081029)[speciation=0]:0.061775,LOC_Os05g47500.1[species=Oryza_sativa]:0.161089)[speciation=1]:0.008388,
Sb09g027330.1[species=Sorghum_bicolor]:0.174412)[speciation=1]:0.105485)[speciation=0]:0.045658)[speciation=0]:0.057427,(Bradi3g10680.1[species=Brachypodium_dista
chyon]:0.342962,((Solyc03g122050.1.1[species=Solanum_lycopersicum]:3.033803,(((IMGA_Medtr2g018510.1[species=Medicago_truncatula]:0.094538,IMGA_Medtr2g018460.1[spe
cies=Medicago_truncatula]:0.194542)[speciation=0]:0.120414,LOC_Os02g21750.1[species=Oryza_sativa]:0.316139)[speciation=1]:0.131131,POPTR_0228s00200.1[species=Popu
lus_trichocarpa]:0.687863)[speciation=0]:0.170176)[speciation=0]:0.258036,(((((((((((((IMGA_Medtr1g086100.1[species=Medicago_truncatula]:0.016825,IMGA_Medtr1g0861
40.1[species=Medicago_truncatula]:0.013388)[speciation=0]:0.279235,GRMZM2G153961_P01[species=Zea_mays]:0.604565)[speciation=1]:0.200312,(Glyma18g24290.1[species=G
lycine_max]:0.187671,IMGA_Medtr6g009080.1[species=Medicago_truncatula]:0.139841)[speciation=1]:0.560529)[speciation=0]:0.113270,LOC_Os05g04610.1[species=Oryza_sat
iva]:0.259335)[speciation=0]:0.230231,IMGA_Medtr2g018480.1[species=Medicago_truncatula]:0.227400)[speciation=0]:0.186191,(evm.model.supercontig_102.37[species=Car
ica_papaya]:0.588904,Solyc03g122070.1.1[species=Solanum_lycopersicum]:0.652276)[speciation=1]:0.068207)[speciation=0]:0.088183,X_29822.m003425[species=Ricinus_com
munis]:0.542472)[speciation=0]:0.212081,IMGA_Medtr6g009070.1[species=Medicago_truncatula]:0.308856)[speciation=0]:0.028462,Glyma12g16410.1[species=Glycine_max]:0.
378669)[speciation=0]:0.084312,Glyma18g01610.1[species=Glycine_max]:0.318807)[speciation=0]:0.067563,GRMZM2G146034_P01[species=Zea_mays]:0.195086)[speciation=0]:0
.154248,(IMGA_Medtr8g066710.1[species=Medicago_truncatula]:0.484493,(GRMZM2G125424_P03[species=Zea_mays]:0.062073,GRMZM2G085236_P03[species=Zea_mays]:0.242416)[sp
eciation=0]:0.333319)[speciation=1]:0.117037)[speciation=0]:0.162561,Solyc05g013890.1.1[species=Solanum_lycopersicum]:0.322329)[speciation=0]:0.204978)[speciation
=0]:0.169881)[speciation=0]:0.072874)[speciation=0]:0.000993)[speciation=0]:0.021118)[speciation=0]:0.017847)[speciation=1]:0.081453)';

#$gene_tree_newick = '(a, (b, ((d,e), (f,g))))';

$gene_tree_newick =~ s/\s+//g;
$gene_tree_newick =~ s/\[speciation=\d\]//g;
print STDERR $gene_tree_newick, "\n";

my $gene_tree = CXGN::Phylo::Parse_newick->new($gene_tree_newick, 0)->parse();
if (!$gene_tree) {
  die"gene tree. Parse_newick->parse() failed to return a tree object. Newick string: "
    . $gene_tree_newick."\n";
}

my $species_newick = "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
my $sparser = CXGN::Phylo::Parse_newick->new($species_newick, 0);
my $species_tree = $sparser->parse();

#my $species_tree = undef;
#my @gtrees = ();
foreach (1..$N) {

 my $orthologger_obj = Orthologger->new({'gene_tree_newick' => $gene_tree_newick, 'species_tree' => $species_tree, 'reroot_method' => 'mindl', 'query_species' => 'castorbean', 'query_id_regex' => undef});

    my $orthologger_outstring = $orthologger_obj->ortholog_result_string();
$orthologger_obj->decircularize();
# find_cycle($orthologger_obj);
# store_orthologger_out($orthologger_outstring, \%idpair_orthocount);
print STDERR "$_\n";
  # my $gtcopy = $gene_tree->get_root()->copy_subtree();
  # #push @gtrees, $gene_tree;
  # #find_cycle($gtcopy);
  # #print STDERR "decircularizing gtcopy.\n";
  # $gtcopy->decircularize();
  # print STDERR "find cycle report: \n";
  # find_cycle($gtcopy);
  # print STDERR "done with tree number: $_.\n";
}
print STDERR "done with all.\n";
