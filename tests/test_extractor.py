import pytest

from hybran import extractor


@pytest.mark.parametrize("description", [
    'clean',
    'interspersed_text',
    'none',
])
def test_repossess(description):
    inputs = {
        'clean': '[topology=circular] [location=chromosome]',
        'interspersed_text': 'random stuff [  topology = circular  ]more random stuff [location=chromosome]',
        'none': None,
    }
    expected = {
        'clean': {
            'topology':'circular',
            'location':'chromosome',
        },
        'interspersed_text': {
            'topology':'circular',
            'location':'chromosome',
        },
        'none': {},
    }

    assert extractor.repossess(inputs[description]) == expected[description]
