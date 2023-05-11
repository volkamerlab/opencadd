"""
test for the module `handle_tautomers`

Partially derived from MolVS's tests and RDKIT MolStandardize tutorial:
https://github.com/mcs07/MolVS/blob/master/tests/test_tautomer.py
"""
import pytest

from rdkit import Chem

from opencadd.compounds.standardization.handle_tautomers import enumerate_tautomer, canonicalize_tautomer


def test_1_3_keto_enol_enumeration():
    """Enumerate 1,3 keto/enol tautomer."""
    assert enumerate_tautomer("C1(=CCCCC1)O") == {"OC1=CCCCC1", "O=C1CCCCC1"}
    assert enumerate_tautomer("C1(CCCCC1)=O") == {"OC1=CCCCC1", "O=C1CCCCC1"}


def test_acetophenone_keto_enol_enumeration():
    """Enumerate acetophenone keto/enol tautomer."""
    assert enumerate_tautomer("C(=C)(O)C1=CC=CC=C1") == {
        "C=C(O)c1ccccc1",
        "CC(=O)c1ccccc1",
    }
    assert enumerate_tautomer("CC(C)=O") == {"CC(C)=O", "C=C(C)O"}


def test_1_5_keto_enol_enumeration3():
    """1,5 keto/enol tautomer"""
    assert enumerate_tautomer("C1(=CC=CCC1)O") == {
        "O=C1C=CCCC1",
        "OC1=CCC=CC1",
        "OC1=CC=CCC1",
        "O=C1CC=CCC1",
        "OC1=CCCC=C1",
    }


def test_aliphatic_imine_enumeration():
    """aliphatic imine tautomer"""
    assert enumerate_tautomer("C1(CCCCC1)=N") == {"N=C1CCCCC1", "NC1=CCCCC1"}
    assert enumerate_tautomer("C1(=CCCCC1)N") == {"N=C1CCCCC1", "NC1=CCCCC1"}


def test_special_imine_enumeration():
    """special imine tautomer"""
    assert enumerate_tautomer("C1(C=CC=CN1)=CC") == {
        "CC=C1C=CC=CN1",
        "CCc1ccccn1",
        "CC=C1C=CCC=N1",
    }
    assert enumerate_tautomer("C1(=NC=CC=C1)CC") == {
        "CC=C1C=CC=CN1",
        "CCc1ccccn1",
        "CC=C1C=CCC=N1",
    }


def test_1_3_aromatic_heteroatom_enumeration():
    """1,3 aromatic heteroatom H shift"""
    assert enumerate_tautomer("O=c1cccc[nH]1") == {"Oc1ccccn1", "O=c1cccc[nH]1"}
    assert enumerate_tautomer("Oc1ccccn1") == {"Oc1ccccn1", "O=c1cccc[nH]1"}
    assert enumerate_tautomer("Oc1ncc[nH]1") == {"Oc1ncc[nH]1", "O=c1[nH]cc[nH]1"}


def test_1_3_heteroatom_enumeration():
    """1,3 heteroatom H shift"""
    assert enumerate_tautomer("OC(C)=NC") == {"CN=C(C)O", "CNC(C)=O", "C=C(O)NC"}
    assert enumerate_tautomer("CNC(C)=O") == {"CN=C(C)O", "CNC(C)=O", "C=C(O)NC"}
    assert enumerate_tautomer("S=C(N)N") == {"N=C(N)S", "NC(N)=S"}
    assert enumerate_tautomer("SC(N)=N") == {"N=C(N)S", "NC(N)=S"}
    assert enumerate_tautomer("N=c1[nH]ccn(C)1") == {"Cn1ccnc1N", "Cn1cc[nH]c1=N"}
    assert enumerate_tautomer("CN=c1[nH]cncc1") == {
        "CN=c1ccnc[nH]1",
        "CNc1ccncn1",
        "CN=c1cc[nH]cn1",
    }


def test_1_5_aromatic_heteroatom_enumeration():
    """1,5 aromatic heteroatom H shift"""
    assert enumerate_tautomer("Oc1cccc2ccncc12") == {
        "O=c1cccc2cc[nH]cc1-2",
        "Oc1cccc2ccncc12",
    }
    assert enumerate_tautomer("O=c1cccc2cc[nH]cc1-2") == {
        "O=c1cccc2cc[nH]cc1-2",
        "Oc1cccc2ccncc12",
    }
    assert enumerate_tautomer("Cc1n[nH]c2ncnn12") == {
        "C=C1NNc2ncnn21",
        "Cc1n[nH]c2ncnn12",
        "Cc1nnc2[nH]cnn12",
        "C=C1NN=C2N=CNN12",
        "Cc1nnc2nc[nH]n12",
        "C=C1NN=C2NC=NN12",
    }
    assert enumerate_tautomer("Cc1nnc2nc[nH]n12") == {
        "C=C1NNc2ncnn21",
        "Cc1n[nH]c2ncnn12",
        "Cc1nnc2[nH]cnn12",
        "C=C1NN=C2N=CNN12",
        "Cc1nnc2nc[nH]n12",
        "C=C1NN=C2NC=NN12",
    }
    assert enumerate_tautomer("Oc1ccncc1") == {"Oc1ccncc1", "O=c1cc[nH]cc1"}
    assert enumerate_tautomer("Oc1c(cccc3)c3nc2ccncc12") == {
        "Oc1c2ccccc2nc2ccncc12",
        "O=c1c2ccccc2[nH]c2ccncc12",
        "O=c1c2c[nH]ccc-2nc2ccccc12",
    }
    assert enumerate_tautomer("C2(=C1C(=NC=N1)[NH]C(=N2)N)O") == {
        "N=c1[nH]c2ncnc-2c(O)[nH]1",
        "Nc1nc(O)c2ncnc-2[nH]1",
        "N=c1nc(O)c2nc[nH]c2[nH]1",
        "Nc1nc2ncnc-2c(O)[nH]1",
        "N=c1nc2nc[nH]c2c(O)[nH]1",
        "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
        "N=c1nc(O)c2[nH]cnc2[nH]1",
        "N=c1[nH]c(=O)c2[nH]cnc2[nH]1",
        "Nc1nc(=O)c2nc[nH]c2[nH]1",
        "Nc1nc(O)c2nc[nH]c2n1",
        "Nc1nc(=O)c2[nH]cnc2[nH]1",
        "N=c1nc2[nH]cnc2c(O)[nH]1",
        "Nc1nc2[nH]cnc2c(=O)[nH]1",
        "Nc1nc2nc[nH]c2c(=O)[nH]1",
        "Nc1nc(O)c2[nH]cnc2n1",
    }
    assert enumerate_tautomer("C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O") == {
        "N=c1[nH]c2ncnc-2c(O)[nH]1",
        "Nc1nc(O)c2ncnc-2[nH]1",
        "N=c1nc(O)c2nc[nH]c2[nH]1",
        "Nc1nc2ncnc-2c(O)[nH]1",
        "N=c1nc2nc[nH]c2c(O)[nH]1",
        "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
        "N=c1nc(O)c2[nH]cnc2[nH]1",
        "N=c1[nH]c(=O)c2[nH]cnc2[nH]1",
        "Nc1nc(=O)c2nc[nH]c2[nH]1",
        "Nc1nc(O)c2nc[nH]c2n1",
        "Nc1nc(=O)c2[nH]cnc2[nH]1",
        "N=c1nc2[nH]cnc2c(O)[nH]1",
        "Nc1nc2[nH]cnc2c(=O)[nH]1",
        "Nc1nc2nc[nH]c2c(=O)[nH]1",
        "Nc1nc(O)c2[nH]cnc2n1",
    }
    assert enumerate_tautomer("Oc1n(C)ncc1") == {
        "Cn1nccc1O",
        "CN1N=CCC1=O",
        "Cn1[nH]ccc1=O",
    }
    assert enumerate_tautomer("O=c1nc2[nH]ccn2cc1") == {
        "O=c1ccn2cc[nH]c2n1",
        "Oc1ccn2ccnc2n1",
        "O=c1ccn2ccnc2[nH]1",
    }
    assert enumerate_tautomer("N=c1nc[nH]cc1") == {
        "N=c1cc[nH]cn1",
        "N=c1ccnc[nH]1",
        "Nc1ccncn1",
    }
    assert enumerate_tautomer("N=c(c1)ccn2cc[nH]c12") == {
        "N=c1ccn2cc[nH]c2c1",
        "Nc1ccn2ccnc2c1",
    }
    assert enumerate_tautomer("CN=c1nc[nH]cc1") == {
        "CN=c1ccnc[nH]1",
        "CNc1ccncn1",
        "CN=c1cc[nH]cn1",
    }


def test_1_3_1_5_aromatic_heteroatom_enumeration():
    """1,3 and 1,5 aromatic heteroatom H shift"""
    assert enumerate_tautomer("Oc1ncncc1") == {
        "Oc1ccncn1",
        "O=c1ccnc[nH]1",
        "O=c1cc[nH]cn1",
    }


def test_1_7_aromatic_heteroatom_enumeration():
    """1,7 aromatic heteroatom H shift"""
    assert enumerate_tautomer("c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1") == {
        "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
        "c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2",
        "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2",
    }
    assert enumerate_tautomer("c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2") == {
        "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
        "c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2",
        "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2",
    }


def test_1_9_aromatic_heteroatom_enumeration():
    """1,9 aromatic heteroatom H shift"""
    assert enumerate_tautomer("CNc1ccnc2ncnn21") == {
        "CN=c1cc[nH]c2ncnn12",
        "CN=c1ccnc2nc[nH]n12",
        "CN=c1ccnc2[nH]cnn12",
        "CNc1ccnc2ncnn12",
    }
    assert enumerate_tautomer("CN=c1ccnc2nc[nH]n21") == {
        "CN=c1ccnc2nc[nH]n12",
        "CN=c1cc[nH]c2ncnn12",
        "CN=c1ccnc2[nH]cnn12",
        "CNc1ccnc2ncnn12",
    }


def test_1_11_aromatic_heteroatom_enumeration():
    """1,11 aromatic heteroatom H shift"""
    assert enumerate_tautomer("Nc1ccc(C=C2C=CC(=O)C=C2)cc1") == {
        "Nc1ccc(C=C2C=CC(=O)C=C2)cc1",
        "N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1",
        "N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
        "N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1",
    }
    assert enumerate_tautomer("N=C1C=CC(=Cc2ccc(O)cc2)C=C1") == {
        "Nc1ccc(C=C2C=CC(=O)C=C2)cc1",
        "N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1",
        "N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
        "N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1",
    }


def test_heterocyclic_enumeration():
    """heterocyclic tautomer"""
    assert enumerate_tautomer("n1ccc2ccc[nH]c12") == {
        "c1c[nH]c2nccc-2c1",
        "c1cnc2[nH]ccc2c1",
    }
    assert enumerate_tautomer("c1cc(=O)[nH]c2nccn12") == {
        "O=c1ccn2cc[nH]c2n1",
        "Oc1ccn2ccnc2n1",
        "O=c1ccn2ccnc2[nH]1",
    }
    assert enumerate_tautomer("c1cnc2c[nH]ccc12") == {
        "c1cc2cc[nH]c2cn1",
        "c1cc2cc[nH]cc-2n1",
    }
    assert enumerate_tautomer("n1ccc2c[nH]ccc12") == {
        "c1cc2[nH]ccc2cn1",
        "c1cc2c[nH]ccc-2n1",
    }
    assert enumerate_tautomer("c1cnc2ccc[nH]c12") == {
        "c1c[nH]c2ccnc-2c1",
        "c1cnc2cc[nH]c2c1",
    }


def test_furanone_enumeration():
    """furanone tautomer"""
    assert enumerate_tautomer("C1=CC=C(O1)O") == {"Oc1ccco1", "O=C1CC=CO1"}
    assert enumerate_tautomer("O=C1CC=CO1") == {"Oc1ccco1", "O=C1CC=CO1"}


def test_keten_ynol_enumeration():
    """keten/ynol tautomer"""
    assert enumerate_tautomer("CC=C=O") == {"CC=C=O", "CC#CO"}
    assert enumerate_tautomer("CC#CO") == {"CC=C=O", "CC#CO"}


def test_ionic_nitro_aci_nitro_enumeration():
    """ionic nitro/aci-nitro tautomer"""
    assert enumerate_tautomer("C([N+](=O)[O-])C") == {
        "CC[N+](=O)[O-]",
        "CC=[N+]([O-])O",
    }
    assert enumerate_tautomer("C(=[N+](O)[O-])C") == {
        "CC[N+](=O)[O-]",
        "CC=[N+]([O-])O",
    }


def test_oxim_nitroso_enumeration():
    """oxim nitroso tautomer"""
    assert enumerate_tautomer("CC(C)=NO") == {"CC(C)N=O", "CC(C)=NO", "C=C(C)NO"}
    assert enumerate_tautomer("CC(C)N=O") == {"CC(C)N=O", "CC(C)=NO", "C=C(C)NO"}
    assert enumerate_tautomer("O=Nc1ccc(O)cc1") == {
        "O=NC1C=CC(=O)C=C1",
        "O=C1C=CC(=NO)C=C1",
        "O=Nc1ccc(O)cc1",
    }
    assert enumerate_tautomer("O=C1C=CC(=NO)C=C1") == {
        "O=NC1C=CC(=O)C=C1",
        "O=C1C=CC(=NO)C=C1",
        "O=Nc1ccc(O)cc1",
    }


def test_cyano_iso_cyanic_acid_enumeration():
    """cyano/iso-cyanic acid tautomer"""
    assert enumerate_tautomer("C(#N)O") == {"N#CO", "N=C=O"}
    assert enumerate_tautomer("C(=N)=O") == {"N#CO", "N=C=O"}


def test_isocyanide_enumeration():
    """isocyanide tautomer"""
    assert enumerate_tautomer("C#N") == {"[C-]#[NH+]", "C#N"}
    assert enumerate_tautomer("[C-]#[NH+]") == {"[C-]#[NH+]", "C#N"}


def test_phosphonic_acid_enumeration():
    """phosphonic acid tautomer"""
    assert enumerate_tautomer("[PH](=O)(O)(O)") == {"OP(O)O", "O=[PH](O)O"}
    assert enumerate_tautomer("P(O)(O)O") == {"OP(O)O", "O=[PH](O)O"}


def test_mobile_double_stereochemistry_enumeration():
    """Remove stereochemistry from mobile double bonds"""
    assert enumerate_tautomer("c1(ccccc1)/C=C(/O)\\C") == {
        "C=C(O)Cc1ccccc1",
        "CC(O)=Cc1ccccc1",
        "CC(=O)Cc1ccccc1",
    }
    assert enumerate_tautomer("C/C=C/C(C)=O") == {
        "C=C(O)C=CC",
        "C=CCC(=C)O",
        "CC=CC(C)=O",
        "C=CCC(C)=O",
        "C=CC=C(C)O",
    }
    assert enumerate_tautomer("C/C=C\\C(C)=O") == {
        "C=C(O)C=CC",
        "C=CCC(=C)O",
        "CC=CC(C)=O",
        "C=CCC(C)=O",
        "C=CC=C(C)O",
    }


def test_gaunine_enumeration():
    """Gaunine tautomers"""
    assert enumerate_tautomer("N1C(N)=NC=2N=CNC2C1=O") == {
        "N=c1[nH]c(=O)c2[nH]cnc2[nH]1",
        "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
        "N=c1[nH]c2ncnc-2c(O)[nH]1",
        "N=c1nc(O)c2[nH]cnc2[nH]1",
        "N=c1nc(O)c2nc[nH]c2[nH]1",
        "N=c1nc2[nH]cnc2c(O)[nH]1",
        "N=c1nc2nc[nH]c2c(O)[nH]1",
        "Nc1nc(=O)c2[nH]cnc2[nH]1",
        "Nc1nc(=O)c2nc[nH]c2[nH]1",
        "Nc1nc(O)c2[nH]cnc2n1",
        "Nc1nc(O)c2nc[nH]c2n1",
        "Nc1nc(O)c2ncnc-2[nH]1",
        "Nc1nc2[nH]cnc2c(=O)[nH]1",
        "Nc1nc2nc[nH]c2c(=O)[nH]1",
        "Nc1nc2ncnc-2c(O)[nH]1",
    }


def test_many_enumeration():
    """Test a structure with hundreds of tautomers."""
    assert len(enumerate_tautomer("[H][C](CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O")) == 375


def test_1_3_keto_enol_canonicalization():
    """1,3 keto/enol tautomer"""
    assert canonicalize_tautomer("C1(=CCCCC1)O") == "O=C1CCCCC1"
    assert canonicalize_tautomer("C1(CCCCC1)=O") == "O=C1CCCCC1"


def test_acetophenone_keto_enol_canonicalization():
    """Acetophenone keto/enol tautomer"""
    assert canonicalize_tautomer("C(=C)(O)C1=CC=CC=C1") == "CC(=O)c1ccccc1"


def test_acetone_keto_enol_canonicalization():
    """Acetone keto/enol tautomer"""
    assert canonicalize_tautomer("CC(C)=O") == "CC(C)=O"


def test_keto_enol_canonicalization():
    """keto/enol tautomer"""
    assert canonicalize_tautomer("OC(C)=C(C)C") == "CC(=O)C(C)C"


def test_phenylpropanone_keto_enol_canonicalization():
    """1-phenyl-2-propanone enol/keto"""
    assert canonicalize_tautomer("c1(ccccc1)CC(=O)C") == "CC(=O)Cc1ccccc1"


def test_1_5_keto_enol_canonicalization():
    """1,5 keto/enol tautomer"""
    assert (
        canonicalize_tautomer("Oc1nccc2cc[nH]c(=N)c12") == "N=c1[nH]ccc2cc[nH]c(=O)c12"
    )
    assert canonicalize_tautomer("C1(C=CCCC1)=O") == "O=C1C=CCCC1"
    assert canonicalize_tautomer("C1(=CC=CCC1)O") == "O=C1C=CCCC1"


def test_aliphatic_imine_canonicalization():
    """aliphatic imine tautomer"""
    assert canonicalize_tautomer("C1(CCCCC1)=N") == "N=C1CCCCC1"
    assert canonicalize_tautomer("C1(=CCCCC1)N") == "N=C1CCCCC1"


def test_special_imine_canonicalization():
    """special imine tautomer"""
    assert canonicalize_tautomer("C1(C=CC=CN1)=CC") == "CCc1ccccn1"
    assert canonicalize_tautomer("C1(=NC=CC=C1)CC") == "CCc1ccccn1"


def test_1_3_aromatic_heteroatom_canonicalization():
    """1,3 aromatic heteroatom H shift"""
    assert canonicalize_tautomer("O=c1cccc[nH]1") == "O=c1cccc[nH]1"
    assert canonicalize_tautomer("Oc1ccccn1") == "O=c1cccc[nH]1"
    assert canonicalize_tautomer("Oc1ncc[nH]1") == "O=c1[nH]cc[nH]1"


def test_1_3_heteroatom_canonicalization():
    """1,3 heteroatom H shift"""
    assert canonicalize_tautomer("OC(C)=NC") == "CNC(C)=O"
    assert canonicalize_tautomer("CNC(C)=O") == "CNC(C)=O"
    assert canonicalize_tautomer("S=C(N)N") == "NC(N)=S"
    assert canonicalize_tautomer("SC(N)=N") == "NC(N)=S"
    assert canonicalize_tautomer("N=c1[nH]ccn(C)1") == "Cn1cc[nH]c1=N"
    assert canonicalize_tautomer("CN=c1[nH]cncc1") == "CN=c1cc[nH]cn1"


def test_1_5_aromatic_heteroatom_canonicalization():
    """1,5 aromatic heteroatom H shift"""
    assert canonicalize_tautomer("Oc1cccc2ccncc12") == "Oc1cccc2ccncc12"
    assert canonicalize_tautomer("O=c1cccc2cc[nH]cc1-2") == "Oc1cccc2ccncc12"
    assert canonicalize_tautomer("Cc1n[nH]c2ncnn12") == "Cc1n[nH]c2ncnn12"
    assert canonicalize_tautomer("Cc1nnc2nc[nH]n12") == "Cc1n[nH]c2ncnn12"
    assert canonicalize_tautomer("Oc1ccncc1") == "O=c1cc[nH]cc1"
    assert (
        canonicalize_tautomer("Oc1c(cccc3)c3nc2ccncc12") == "O=c1c2ccccc2[nH]c2ccncc12"
    )
    assert (
        canonicalize_tautomer("C2(=C1C(=NC=N1)[NH]C(=N2)N)O")
        == "N=c1[nH]c(=O)c2[nH]cnc2[nH]1"
    )
    assert (
        canonicalize_tautomer("C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O")
        == "N=c1[nH]c(=O)c2[nH]cnc2[nH]1"
    )
    assert canonicalize_tautomer("Oc1n(C)ncc1") == "Cn1[nH]ccc1=O"
    assert canonicalize_tautomer("O=c1nc2[nH]ccn2cc1") == "O=c1ccn2cc[nH]c2n1"
    assert canonicalize_tautomer("N=c1nc[nH]cc1") == "N=c1cc[nH]cn1"
    assert canonicalize_tautomer("N=c(c1)ccn2cc[nH]c12") == "N=c1ccn2cc[nH]c2c1"
    assert canonicalize_tautomer("CN=c1nc[nH]cc1") == "CN=c1cc[nH]cn1"


def test_1_3_1_5_aromatic_heteroatom_canonicalization():
    """1,3 and 1,5 aromatic heteroatom H shift"""
    assert canonicalize_tautomer("Oc1ncncc1") == "O=c1cc[nH]cn1"


def test_1_7_aromatic_heteroatom_canonicalization():
    """1,7 aromatic heteroatom H shift"""
    assert (
        canonicalize_tautomer("c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1")
        == "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1"
    )
    assert (
        canonicalize_tautomer("c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2")
        == "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1"
    )


def test_1_9_aromatic_heteroatom_canonicalization():
    """1,9 aromatic heteroatom H shift"""
    assert canonicalize_tautomer("CNc1ccnc2ncnn21") == "CN=c1cc[nH]c2ncnn12"
    assert canonicalize_tautomer("CN=c1ccnc2nc[nH]n21") == "CN=c1cc[nH]c2ncnn12"


def test_1_11_aromatic_heteroatom_canonicalization():
    """1,11 aromatic heteroatom H shift"""
    assert (
        canonicalize_tautomer("Nc1ccc(C=C2C=CC(=O)C=C2)cc1")
        == "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"
    )
    assert (
        canonicalize_tautomer("N=C1C=CC(=Cc2ccc(O)cc2)C=C1")
        == "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"
    )


def test_heterocyclic_canonicalization():
    """heterocyclic tautomer"""
    assert canonicalize_tautomer("n1ccc2ccc[nH]c12") == "c1cnc2[nH]ccc2c1"
    assert canonicalize_tautomer("c1cc(=O)[nH]c2nccn12") == "O=c1ccn2cc[nH]c2n1"
    assert canonicalize_tautomer("c1cnc2c[nH]ccc12") == "c1cc2cc[nH]c2cn1"
    assert canonicalize_tautomer("n1ccc2c[nH]ccc12") == "c1cc2[nH]ccc2cn1"
    assert canonicalize_tautomer("c1cnc2ccc[nH]c12") == "c1cnc2cc[nH]c2c1"


def test_furanone_canonicalization():
    """furanone tautomer"""
    assert canonicalize_tautomer("C1=CC=C(O1)O") == "Oc1ccco1"
    assert canonicalize_tautomer("O=C1CC=CO1") == "Oc1ccco1"


def test_keten_ynol_canonicalization():
    """keten/ynol tautomer"""
    assert canonicalize_tautomer("CC=C=O") == "CC=C=O"
    assert canonicalize_tautomer("CC#CO") == "CC=C=O"


def test_ionic_nitro_aci_nitro_canonicalization():
    """ionic nitro/aci-nitro tautomer"""
    assert canonicalize_tautomer("C([N+](=O)[O-])C") == "CC[N+](=O)[O-]"
    assert canonicalize_tautomer("C(=[N+](O)[O-])C") == "CC[N+](=O)[O-]"


def test_oxim_nitroso_canonicalization():
    """oxim nitroso tautomer"""
    assert canonicalize_tautomer("CC(C)=NO") == "CC(C)=NO"
    assert canonicalize_tautomer("CC(C)N=O") == "CC(C)=NO"


def test_oxim_nitroso_phenol_canonicalization():
    """oxim/nitroso tautomer via phenol"""
    assert canonicalize_tautomer("O=Nc1ccc(O)cc1") == "O=Nc1ccc(O)cc1"
    assert canonicalize_tautomer("O=C1C=CC(=NO)C=C1") == "O=Nc1ccc(O)cc1"


def test_cyano_iso_cyanic_acid_canonicalization():
    """cyano/iso-cyanic acid tautomer"""
    assert canonicalize_tautomer("C(#N)O") == "N=C=O"
    assert canonicalize_tautomer("C(=N)=O") == "N=C=O"


def test_formamidinesulfinic_acid_canonicalization():
    """formamidinesulfinic acid tautomer"""
    assert canonicalize_tautomer("N=C(N)S(=O)O") == "N=C(N)S(=O)O"


def test_isocyanide_canonicalization():
    """isocyanide tautomer"""
    assert canonicalize_tautomer("C#N") == "C#N"
    assert canonicalize_tautomer("[C-]#[NH+]") == "C#N"


def test_phosphonic_acid_canonicalization():
    """phosphonic acid tautomer"""
    assert canonicalize_tautomer("[PH](=O)(O)(O)") == "O=[PH](O)O"
    assert canonicalize_tautomer("P(O)(O)O") == "O=[PH](O)O"
