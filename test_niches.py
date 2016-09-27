# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 10:49:39 2016

@author: adam
"""
from niches import random_niches, T

def test_random_niches():
    obs = random_niches(10)
    assert isinstance(obs, set)
    assert "mine" in obs
    assert None not in obs
    assert len(obs) <= 10
    for niche in obs:
        assert niche in T
