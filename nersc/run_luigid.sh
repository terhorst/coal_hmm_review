#!/bin/bash
BASE=/project/projectdirs/m2871/terhorst/coal_hmm_review/luigid
mkdir -p $BASE/log
luigid --pidfile $BASE/luigid.pid --state-path $BASE/state
