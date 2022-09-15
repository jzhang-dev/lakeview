import doctest
import lakeview as lv

def test_util():
    doctest.testmod(lv.util)

def test_helpers():
    doctest.testmod(lv.helpers)