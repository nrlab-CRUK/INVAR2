#!/bin/bash

grep 'INVAR_SCORE' GLRT.out | \
sed 's/"//g' | \
sed 's/^\[1\] //' | \
sed 's/iteration:/iteration\n/' | \
sed 's/^ //' | \
sed 's/ /,/g' > /tmp/ivz

head -n 1 /tmp/ivz > /tmp/ivy
sed -n '2~2p' /tmp/ivz >> /tmp/ivy

mv /tmp/ivy $1
rm -f /tmp/iv?

