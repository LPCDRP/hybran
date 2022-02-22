import pytest

from hybran import designator


@pytest.mark.parametrize("qualifiers,expected", [
    (dict(), dict(note=["new_note"])),
    (dict(note=["original_note"]), dict(note=["original_note","new_note"])),
])
def test_append_qualifier(qualifiers, expected):
    designator.append_qualifier(qualifiers,'note','new_note')
    assert qualifiers == expected
