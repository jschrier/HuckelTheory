(* ::Package:: *)

BeginPackage["HuckelTheory`"]

HuckelHamiltonianMatrix::usage = 
	"HuckelHamiltonianMatrix[m] gives the pi-electron Hamiltonian matrix for the Molecule[] m."
     
HuckelMO::usage = 
	"HuckelMO[m] returns an association whose keys are the eigenvalues and whose values are the corresponding eigenvectors (molecular orbitals) of Molecule[] m"

HuckelChargeDensityMatrix::usage = 
	"HuckelChargeDensityMatrix[m] returns the pi-charge density matrix of Molecule[] m"

HuckelMOPlot::usage = 
	"HuckelMOPlot[m] returns an association whose keys are the eigenvalues and whose values are plots of the pi-molecular orbitals on the 2D graph of Molecule[] m"

HuckelBondOrderPlot::usage =
	"HuckelBondOrderPlot[m] returns an image with the nearest-neighbor atom pi-bond orders plotted on the 2D dimensional graph of Molecule[] m"
	
HuckelTotalElectronsPlot::usage = 
	"HuckelTotalElectronsPlot[m] returns an image of Molecule[] m depicting the total number of pi electrons on each atomic site."
	
HuckelNetChargePlot::usage =
	"HuckelNetChargePlot[m] returns an image of Molecule[] m depicting the net pi-charge on each atomic site."

Begin["`Private`"]

Once[If[$VersionNumber<12.,Print["Mathematica version 12 or greater is required."];Abort[];]] 

HuckelHamiltonianMatrix[m_Molecule,t_:-1]:=With[{m2=sortPiAtoms[m]},
With[
{atomTypes=MoleculeProperty["AtomList"]@m2,
piElectrons=effectivePiElectrons@m2,
conjugatedAtomIndices=piSystemIndices@m2,
piBonds=piBondList@m2},
t*(diagonal[atomTypes,conjugatedAtomIndices,piElectrons]+(*diagonal terms*)
offDiagonal[atomTypes,piElectrons,Length[conjugatedAtomIndices],piBonds])
]]

HuckelMO[m_Molecule,t_:-1.]:=Map[Chop,AssociationThread@@Eigensystem@HuckelHamiltonianMatrix[m,t]]//KeySort

HuckelChargeDensityMatrix[m_Molecule,t_:-1.]:=With[
{coeffs=Values@HuckelMO[m,t],
nOcc=effectivePiElectronCount[m]/2, (*todo: add an option for net charge*)
nPiAtoms=Length@piSystemIndices@m},
Table[
Sum[2*coeffs[[i,r]]*coeffs[[i,s]],{i,nOcc}],
{r,nPiAtoms},{s,nPiAtoms}]//Chop]

onSp2Carbon[pattern_,bondType_:"Single"][m_Molecule]:=
With[
{sp2Carbons=AtomList[m,
MoleculePattern[{Atom["C","OrbitalHybridization"->"SP2"]}],
"AtomIndex"],
targetAtoms=Flatten@BondList[m,
MoleculePattern[{pattern,Atom["C","OrbitalHybridization"->"SP2"]},{Bond[{1,2},bondType]}],
"AtomIndex"]},
Complement[targetAtoms,sp2Carbons]]

nonRingSP2Atom[m_Molecule,symbol_String]:=AtomList[
m,
MoleculePattern[{Atom[symbol,"OrbitalHybridization"->"SP2","RingMemberQ"->False]}], 
"AtomIndex"]

piSystemIndices[m_Molecule]:=With[
{sp2Atoms=AtomList[m,
MoleculePattern[{Atom[_,"OrbitalHybridization"->"SP2"]}], (*any sp2 species*)
"AtomIndex"],
halogens=onSp2Carbon["F"|"Br"|"Cl"][m]},
Union[sp2Atoms,halogens]]

sortPiAtoms[m_Molecule]:=With[
{piAtoms=piSystemIndices[m],
allAtoms=MoleculeProperty["AtomIndex"]@m},
MoleculeModify[m, {"RenumberAtoms",Join[piAtoms,Complement[allAtoms,piAtoms]]}]
]

effectivePiElectrons[m_Molecule]:=With[
{default=MoleculeProperty["PiElectronCount"]@m,
conjugatedHalogens=onSp2Carbon["F"|"Br"|"Cl"]@m,
alcoholOrEtherOxygen=onSp2Carbon["O"]@m,
anilineLikeNitrogen=nonRingSP2Atom[m,"N"]
(*need to treat other cases?*)},
ReplacePart[default,Partition[
   Join[conjugatedHalogens,alcoholOrEtherOxygen,anilineLikeNitrogen],
   1]->2]]

effectivePiElectronCount[m_Molecule]:=Total@effectivePiElectrons[m]

piBondList[m_Molecule]:=With[
{bonds=BondList[m],
piAtoms=piSystemIndices[m]},
Select[bonds,ContainsOnly[First[#],piAtoms]&]]

haa[Atom["B"],0]=-1.;
haa[Atom["N"],1]=0.5;
haa[Atom["N"],2]=1.5;
haa[Atom["O"],1]=1.0;
haa[Atom["O"],2]=2.0;
haa[Atom["C"],1]=0.0;
haa[Atom["C"],2]=0.0; (*acetylene should be fine also*)
haa[Atom["F"],2]=3.0;
haa[Atom["Cl"],2]=2.0;
haa[Atom["Br"],2]=1.5;
hab[{Atom["B"],Atom["C"]},_,_]=0.7;
hab[{Atom["B"],Atom["N"]},_,_]=0.8;
hab[{Atom["C"],Atom["N"]} ,_,1]=1.0;
hab[{Atom["C"],Atom["N"]} ,_,2]=0.8;
hab[{Atom["C"],Atom["C"]},_,_]=1.0;
hab[{Atom["C"],Atom["O"]},"Double",1]=1.0; (*carbonyl*)
hab[{Atom["C"],Atom["O"]},_,2]=0.8;
hab[{Atom["C"],Atom["Cl"]},_,_]=0.4;
hab[{Atom["Br"],Atom["C"]},_,_]=0.3;
hab[{Atom["C"],Atom["F"]},_,_]=0.7;

bondType[Bond[pair_List,type_String],atoms_List,piElectrons_List]:=With[
{symbols=Sort@atoms[[pair]],
atom1=atoms[[pair[[1]]]]},
pair->hab[symbols,type,
If[atom1===Atom["C"], (*get the non-carbon number of pi electrons*)
piElectrons[[pair[[2]]]],
piElectrons[[pair[[1]] ]] 
]
]
]

offDiagonal[atomTypes_,piElectrons_,conjugatedAtomCount_,piBonds_]:=With[
{half=SparseArray[
bondType[#,atomTypes,piElectrons]&/@piBonds,{conjugatedAtomCount,conjugatedAtomCount} (*force it to be a square matrix of this size; otherwise it breaks with linear molecules*)
]
},
half+Transpose[half]]

diagonal[atomTypes_,conjugatedAtomIndices_,piElectrons_]:=With[
{rules=MapIndexed[{First@#2,First@#2}->haa@@#1&,Transpose[{atomTypes[[conjugatedAtomIndices]],piElectrons[[conjugatedAtomIndices]]}]]
}, (*warning...this assumes that the conjugated atoms all come first before non-conjugated...in any case this will get screwy with display too*)
SparseArray[rules]]

huckelAdjacencyMatrix[m_Molecule]:=With[
{atomTypes=MoleculeProperty["AtomList"]@m,
piElectrons=effectivePiElectrons@m,
conjugatedAtomCount=Length@piSystemIndices@m,
piBonds=piBondList@m},
(*kind of overkill, but easy*)
Unitize@offDiagonal[atomTypes,piElectrons,conjugatedAtomCount,piBonds]]

HuckelMOPlot[m_Molecule,t_:-1.]:=With[
{MOs=N@HuckelMO[m,t],
nOcc=effectivePiElectronCount[m]/2}, (*assumes all doubly occupied MOs*)
AssociationThread[
Keys[MOs]-> 
MapThread[
moPlot[#1][m,#2]&,
{ConstantArray[occupiedMO,nOcc]~Join~ConstantArray[unoccupiedMO,(Length[MOs]-nOcc)],
Values@MOs}]
]]

HuckelTotalElectronsPlot[m_Molecule,t_:-1.]:=With[
{nElectronsHuckel=Chop@Diagonal@HuckelChargeDensityMatrix[m,t]},
moPlot[netCharge][m,nElectronsHuckel]]

HuckelNetChargePlot[m_Molecule,t_:-1.]:=With[
{nElectronsHuckel=Chop@Diagonal@HuckelChargeDensityMatrix[m,t],
nElectronsInitial=effectivePiElectrons[m]},
moPlot[netCharge][m,(*just take the subset of pi orbitals*)
nElectronsInitial[[1;;Length[nElectronsHuckel]]]-nElectronsHuckel] 
]

HuckelBondOrderPlot[m_Molecule,t_:-1.]:=With[
{chargeDensityMatrix=Chop@HuckelChargeDensityMatrix[m,t],
piAdjacencyMatrix=huckelAdjacencyMatrix[m]},
boPlot[bondOrder][m,chargeDensityMatrix*piAdjacencyMatrix]]

moPlot[colorRule_][m_Molecule,coeffs_List]:=With[
{m2=MoleculeModify["RemoveHydrogens"]@sortPiAtoms[m]},
With[
{xyCoords=m2["AtomDiagramCoordinates"][[piSystemIndices[m2]]],
refArea=(1/Max[ Max[#,10^-2]&/@Abs[coeffs]]) }, (*handle cases where charge is zero?*)
Show[
MoleculePlot[m2,PlotTheme->"HeavyAtom"],
Graphics@MapThread[
colorRule[#1,#2,refArea]&,(*turn off normalization scaling*)
{coeffs,xyCoords}]]
]]

boPlot[colorRule_][m_Molecule,coeffs_?MatrixQ]:=With[
{m2=MoleculeModify["RemoveHydrogens"]@sortPiAtoms[m]},
With[
{xyCoords=m2["AtomDiagramCoordinates"][[piSystemIndices[m2]]],
nAtoms=Length[coeffs]},
Show[
MoleculePlot[m2,PlotTheme->"HeavyAtom"],
Graphics@Flatten@Table[
colorRule[coeffs[[r,s]],Midpoint[xyCoords[[{r,s}]]]],
{r,1,nAtoms-1},{s,r+1,nAtoms}]
]
]]

(*I thought I would want independent control over display characteristics like radius, but that's over-rated.  Could probably be refactored*)
occupiedMO[coeff_?NumericQ,xy_List,scaleFactor_:1,arbitraryRadius_:0.5]:=Tooltip[{If[coeff<0,Red,Blue],Opacity[0.3],Disk[xy,Abs[coeff*scaleFactor]*arbitraryRadius]},
coeff]

unoccupiedMO[coeff_?NumericQ,xy_List,scaleFactor_:1.,arbitraryRadius_:0.5]:=Tooltip[{If[coeff<0,Yellow,Green],Opacity[0.4],Disk[xy,Abs[coeff*scaleFactor]*arbitraryRadius]},
coeff]

netCharge[coeff_?NumericQ,xy_List,scaleFactor_:1., arbitraryRadius_:0.5]:=Tooltip[
{If[coeff<0,Red,Black],Opacity[0.3],Disk[xy,Abs[coeff*scaleFactor]*arbitraryRadius]},
coeff]

bondOrder[coeff_?NumericQ,xy_List,scaleFactor_:1.,arbitraryRadius_:0.5]:=
Tooltip[
{Orange,Opacity[0.5],Disk[xy,Abs[coeff*scaleFactor]*arbitraryRadius]},
coeff]

End[ ]

EndPackage[ ]






