import pytest

from hybran import designator


@pytest.mark.parametrize("qualifiers,expected", [
    (dict(), dict(note=["new_note"])),
    (dict(note=["original_note"]), dict(note=["original_note","new_note"])),
])
def test_append_qualifier(qualifiers, expected):
    designator.append_qualifier(qualifiers,'note','new_note')
    assert qualifiers == expected


@pytest.mark.parametrize("name,expected", [
    ('Rv1234::ORF1234', True),
    ('ORF1234', True),
])
def test_has_unannotated_component(name, expected):
    designator.generic_orf_prefix = ['ORF']

    assert designator.has_unannotated_component(name) == expected
