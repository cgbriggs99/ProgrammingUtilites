#!/usr/bin/python3

import os
import sys
import pytest

__tests = ["test_roots", "test_utils"]

def test_cfuncs() :
    for test in __tests :
        os.system("./ctests/" + test)

def test_memory() :
    for test in __tests :
        assert(os.system("valgrind --error-exitcode=1 ./ctests/" + test) != 1)


    
