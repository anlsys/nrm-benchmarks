#!/bin/sh
set -ex

# basic setup
mkdir -p m4
aclocal -I m4
autoreconf -fi
