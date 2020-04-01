import pytest
from .mmligner import MMLignerWrapper
from .base import CommandLineWrapper


@pytest.mark.parametrize(
    "structure1, structure2, structure3",
    [
        "ARN",
        "DBC",
        "EQZGHILKM"
    ],
)
def test_calculate_invalid_inputs(structure1, structure2, structure3):
    mmligner = MMLignerWrapper

    with pytest.raises(TypeError):
        mmligner.calculate([[structure1], [None]])

    with pytest.raises(TypeError):
        mmligner.calculate([None, structure2])

     with pytest.raises(TypeError):
         mmligner.calculate([structure1])

    with pytest.raises(TypeError):
        mmligner.calculate([structure1, structure2, structure3])


@pytest.mark.parametrize(
    "structure1, structure2, result",
    [
        "SOME_STRUCTURE", #TODO: Find good test sample(s)
        "SOME_STRUCTURE",
        {
            "RMSD" = ,
            "Score" = ,
            "metadata" = {}
        },
    ],
)
def test_calculate_results(structure1, structure2, result)
    mmligner = MMLignerWrapper

    calculate_result = mmligner.calculate([structure1, structure2])

    assert pytest.approx(calculate_result["RMSD"], result["RMSD"])
    assert pytest.approx(calculate_result["Score"], result["Score"])
    assert(calculate_result("metadata")==result["metadata"])


def test_deletion_of_temp_files() #Do we need something like that?
