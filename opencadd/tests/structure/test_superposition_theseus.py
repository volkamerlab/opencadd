"""
Test opencadd.structure.superposition.theseus
"""

import os
import pytest
from opencadd.structure.superposition.api import Structure
from opencadd.structure.superposition.engines.theseus import TheseusAligner


def test_theseus_instantiation():
    aligner = TheseusAligner()


@pytest.mark.skipif(
    "GITHUB_ACTIONS" in os.environ, reason="FIXME: GH Actions segfaults with this test"
)
def test_theseus_different_no_superposition():
    different_models = [Structure.from_pdbid(pdb_id) for pdb_id in ["2BBM", "1CFC"]]
    aligner = TheseusAligner(statistics_only=True)
    results = aligner.calculate(different_models)
    assert "scores" in results
    assert "rmsd" in results["scores"]
    assert "coverage" in results["scores"]
    assert "metadata" in results
    assert "least_squares" in results["metadata"]
    assert "maximum_likelihood" in results["metadata"]
    assert "log_marginal_likelihood" in results["metadata"]
    assert "aic" in results["metadata"]
    assert "bic" in results["metadata"]
    assert "omnibus_chi_square" in results["metadata"]
    assert "hierarchical_var_chi_square" in results["metadata"]
    assert "rotational_translational_covar_chi_square" in results["metadata"]
    assert "hierarchical_minimum_var" in results["metadata"]
    assert "hierarchical_minimum_sigma" in results["metadata"]
    assert "skewness" in results["metadata"]
    assert "skewness_z" in results["metadata"]
    assert "kurtosis" in results["metadata"]
    assert "kurtosis_z" in results["metadata"]
    assert "data_pts" in results["metadata"]
    assert "free_params" in results["metadata"]
    assert "d_p" in results["metadata"]
    assert "median_structure" in results["metadata"]
    assert "n_total" in results["metadata"]
    assert "n_structures" in results["metadata"]
    assert "total_rounds" in results["metadata"]


def test_theseus_different():
    different_models = [Structure.from_pdbid(pdb_id) for pdb_id in ["6HG4", "6HG9"]]
    aligner = TheseusAligner()
    results = aligner.calculate(different_models, identical=False)
    assert "superposed" in results
    assert "scores" in results
    assert "rmsd" in results["scores"]
    assert "coverage" in results["scores"]
    assert "metadata" in results
    assert "transformation" in results["metadata"]
    assert "least_squares" in results["metadata"]
    assert "maximum_likelihood" in results["metadata"]
    assert "log_marginal_likelihood" in results["metadata"]
    assert "aic" in results["metadata"]
    assert "bic" in results["metadata"]
    assert "omnibus_chi_square" in results["metadata"]
    assert "hierarchical_var_chi_square" in results["metadata"]
    assert "rotational_translational_covar_chi_square" in results["metadata"]
    assert "hierarchical_minimum_var" in results["metadata"]
    assert "hierarchical_minimum_sigma" in results["metadata"]
    assert "skewness" in results["metadata"]
    assert "skewness_z" in results["metadata"]
    assert "kurtosis" in results["metadata"]
    assert "kurtosis_z" in results["metadata"]
    assert "data_pts" in results["metadata"]
    assert "free_params" in results["metadata"]
    assert "d_p" in results["metadata"]
    assert "median_structure" in results["metadata"]
    assert "n_total" in results["metadata"]
    assert "n_structures" in results["metadata"]
    assert "total_rounds" in results["metadata"]
    assert pytest.approx(results["scores"]["rmsd"], 1.54)


def test_metadata_parsing():
    output = """
                            < BEGIN THESEUS 3.3.0 >
I===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-=I
I                THESEUS: Maximum likelihood multiple superpositioning        I
I=-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===I
    Detected 8 CPUs ...
    Reading 2 pdb files ...
    Successfully read 2 models and/or structures
    Reading multiple sequence alignment ...
    Calculating superposition transformations ...
    Calculating statistics ...
    Calculating likelihood statistics ...
    2 models superimposed in 6.3 ms
  * Classical LS pairwise <RMSD>                    1.54018
  * Least-squares <sigma>                           0.44461
  * Maximum Likelihood <sigma>                      0.27576
  ~ Marginal Log Likelihood                        -1394.25
  ~ AIC                                            -4610.73
  ~ BIC                                            -7817.24
  + Omnibus chi^2                                      5.02 (P:0.00e+00)
  + Hierarchical var (1.14e-01, 1.50e+00) chi^2        1.85 (P:1.00e+00)
  + Rotational, translational, covar chi^2             5.02 (P:0.00e+00)
  + Hierarchical minimum var (sigma)               2.07e-02 (1.44e-01)
  < skewness                                           0.00 (P:1.00e+00)
  < skewness Z-value                                   0.00
  < kurtosis                                          -0.86 (P:8.37e-23)
  < kurtosis Z-value                                   9.83
  * Data pts = 3162,  Free params = 1594,  D/P = 2.0
  * Median structure = #2
  * N(total) = 1054, N(atoms) = 527, N(structures) = 2
  Total rounds = 20
  Converged to a fractional precision of 7.3e-08
I===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-==I
    Transforming coordinates ...
    Writing transformations file ...
    Writing transformed coordinates PDB file ...
    Writing average coordinate file ...
    Done.
I===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-===-==I
                            <  END THESEUS 3.3.0  >
    """

    aligner = TheseusAligner()
    r = aligner._parse_superposition(output)
    assert r["metadata"]["aic"] == -4610.73
