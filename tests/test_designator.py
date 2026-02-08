import pytest

from hybran import designator


@pytest.mark.parametrize("qualifiers,expected", [
    (dict(), dict(note=["new_note"])),
    (dict(note=["original_note"]), dict(note=["original_note","new_note"])),
])
def test_append_qualifier(qualifiers, expected):
    designator.append_qualifier(qualifiers,'note','new_note')
    assert qualifiers == expected


@pytest.mark.parametrize("ids,expected", [
    ({"gene1", "gene2", "gene3"}, 1),
    ({"gene1", "gene2", "HYBRA0001"}, 2),
    ({"HYBRA0004", "gene1", "HYBRA0002", "gene2"}, 5),
    ({"gene1::HYBRA0003", "gene2", "HYBRA0001"}, 4),
])
def test_find_next_increment(ids, expected):
    assert designator.find_next_increment(ids) == expected


@pytest.mark.parametrize("name,expected", [
    ('Rv1234::ORF1234', True),
    ('ORF1234', True),
])
def test_has_unannotated_component(name, expected):
    designator.generic_orf_prefix = ['ORF']

    assert designator.has_unannotated_component(name) == expected
