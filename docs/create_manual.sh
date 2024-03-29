#!/bin/bash

sphinx-build -M latexpdf source build

mv build/latex/waveproshm.pdf manuals/waveproshm.pdf
